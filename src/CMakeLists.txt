add_subdirectory(contrafold)
add_subdirectory(probconsRNA)
add_subdirectory(contralign)
add_subdirectory(linearfold)
add_subdirectory(linearalign)

# for building a statically linked binary
# set(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
# set(BUILD_SHARED_LIBRARIES OFF)
# set(CMAKE_EXE_LINKER_FLAGS "-static")

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

add_executable(
    dafs
    alifold.cpp align.cpp boltzmann_param.c dafs.cpp fa.cpp fold.cpp
    gradient_manager.cpp ip.cpp ipknot.cpp linfold_wrapper.cpp linearalign.cpp needleman_wunsch.cpp nussinov.cpp ribosum.cpp tree.cpp
)

target_compile_definitions(dafs PRIVATE PACKAGE="dafs")
target_compile_definitions(dafs PRIVATE PACKAGE_VERSION="0.0.4")
target_compile_definitions(dafs PRIVATE VERSION="0.0.4")
target_include_directories(dafs PRIVATE "spdlog/include")
target_include_directories(dafs PRIVATE "cxxopts/include")
target_compile_features(dafs PUBLIC cxx_std_17)

find_package(PkgConfig REQUIRED)
pkg_check_modules(VIENNARNA REQUIRED RNAlib2>=2.7)
if (VIENNARNA_FOUND)
  # Only nearest-neighbor energy primitives are used by the linear profile
  # model; its partition/inside/outside algorithms remain in LinFold.
  target_include_directories(linearfold PUBLIC ${VIENNARNA_INCLUDE_DIRS})
  target_link_libraries(linearfold PUBLIC ${VIENNARNA_LDFLAGS} ${CMAKE_DL_LIBS})
  target_include_directories(dafs PRIVATE ${VIENNARNA_INCLUDE_DIRS})
  target_compile_definitions(dafs PRIVATE HAVE_VIENNA18 HAVE_VIENNA20)
  target_link_libraries(dafs PRIVATE ${VIENNARNA_LDFLAGS} ${CMAKE_DL_LIBS})
endif()

if(ENABLE_CPLEX)
  find_package(CPLEX)
  if (CPLEX_FOUND)
    set(DAFS_HAS_IP_BACKEND TRUE)
    target_compile_definitions(dafs PRIVATE WITH_CPLEX)
    target_include_directories(dafs PRIVATE ${CPLEX_INCLUDE_DIR})
    target_link_libraries(dafs PRIVATE ${CPLEX_LIBRARY})
  endif()
elseif(ENABLE_GUROBI)
  find_package(GUROBI)
  if (GUROBI_FOUND)
    set(DAFS_HAS_IP_BACKEND TRUE)
    target_compile_definitions(dafs PRIVATE WITH_GUROBI)
    target_include_directories(dafs PRIVATE ${GUROBI_INCLUDE_DIR})
    target_link_libraries(dafs PRIVATE ${GUROBI_LIBRARY})
  endif()
elseif(ENABLE_SCIP)
  find_package(SCIP)
  if (SCIP_FOUND)
    set(DAFS_HAS_IP_BACKEND TRUE)
    message(STATUS "Found SCIP version ${SCIP_VERSION}")
    target_compile_definitions(dafs PRIVATE WITH_SCIP)
    target_include_directories(dafs PRIVATE ${SCIP_INCLUDE_DIRS})
    target_link_libraries(dafs PRIVATE ${SCIP_LIBRARIES})
  endif()
elseif(ENABLE_HIGHS)
  find_package(HiGHS)
  find_package(ZLIB 1.2.3)
  if (HiGHS_FOUND)
    set(DAFS_HAS_IP_BACKEND TRUE)
    target_compile_definitions(dafs PRIVATE WITH_HIGHS)
    target_include_directories(dafs PRIVATE ${HiGHS_INCLUDE_DIR})
    target_link_libraries(dafs PRIVATE ${HiGHS_LIBRARY} ${ZLIB_LIBRARY})
  endif()
else()
  find_package(GLPK)
  if (GLPK_FOUND)
    set(DAFS_HAS_IP_BACKEND TRUE)
    target_compile_definitions(dafs PRIVATE WITH_GLPK)
    target_include_directories(dafs PRIVATE ${GLPK_INCLUDE_DIR})
    target_link_libraries(dafs PRIVATE ${GLPK_LIBRARY} ${ZLIB_LIBRARY} ${GMP_LIBRARY} ${DL_LIBRARY} ${LTDL_LIBRARY} ${MPFR_LIBRARY})
  endif()
endif()

target_link_libraries(dafs PRIVATE contrafold probcons contralign linearfold linearalign)

if(BUILD_TESTING)
  add_executable(dd_certificate_test ../tests/dd_certificate_test.cpp)
  target_include_directories(dd_certificate_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(dd_certificate_test PRIVATE cxx_std_17)
  add_test(NAME dd_certificate COMMAND dd_certificate_test)
  find_package(Python3 COMPONENTS Interpreter REQUIRED)
  add_test(NAME dd_certificate_pipeline COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_dd_certificates.py $<TARGET_FILE:dafs>)
  add_executable(
      linear_models_fail_fast_test
      ../tests/linear_models_fail_fast.cpp
      linfold_wrapper.cpp linearalign.cpp ribosum.cpp)
  target_include_directories(linear_models_fail_fast_test PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(linear_models_fail_fast_test PRIVATE cxx_std_17)
  target_link_libraries(linear_models_fail_fast_test PRIVATE linearfold linearalign)

  add_test(NAME linear_models_fail_fast COMMAND linear_models_fail_fast_test)

  add_executable(linear_profile_oracle_test
      ../tests/linear_profile_oracle_test.cpp linfold_wrapper.cpp ribosum.cpp)
  target_include_directories(linear_profile_oracle_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(linear_profile_oracle_test PRIVATE cxx_std_17)
  target_link_libraries(linear_profile_oracle_test PRIVATE linearfold)
  add_test(NAME linear_profile_oracle COMMAND linear_profile_oracle_test)

  add_executable(linear_partition_test ../tests/linear_partition_test.cpp)
  target_include_directories(linear_partition_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(linear_partition_test PRIVATE cxx_std_17)
  target_link_libraries(linear_partition_test PRIVATE linearfold)
  add_test(NAME linear_partition_exhaustive COMMAND linear_partition_test)

  add_executable(linear_capacity_test ../tests/linear_capacity_test.cpp)
  target_include_directories(linear_capacity_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(linear_capacity_test PRIVATE cxx_std_17)
  target_link_libraries(linear_capacity_test PRIVATE linearfold)
  add_test(NAME linear_partition_capacity COMMAND linear_capacity_test
      ${PROJECT_SOURCE_DIR}/tests/data/linear_capacity_sequence.fa)

  add_executable(relaxed_bounds_test ../tests/relaxed_bounds_test.cpp)
  target_include_directories(relaxed_bounds_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(relaxed_bounds_test PRIVATE cxx_std_17)
  add_test(NAME relaxed_bounds COMMAND relaxed_bounds_test)

  add_executable(
      gradient_manager_test
      ../tests/gradient_manager_test.cpp gradient_manager.cpp)
  target_include_directories(gradient_manager_test PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(gradient_manager_test PRIVATE cxx_std_17)
  add_test(NAME gradient_manager COMMAND gradient_manager_test)

  add_executable(dd_recovery_test ../tests/dd_recovery_test.cpp)
  target_include_directories(dd_recovery_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(dd_recovery_test PRIVATE cxx_std_17)
  add_test(NAME dd_recovery COMMAND dd_recovery_test)
  add_test(NAME dd_recovery_pipeline COMMAND ${Python3_EXECUTABLE}
      ${PROJECT_SOURCE_DIR}/tests/test_dd_recovery.py $<TARGET_FILE:dafs>)

  add_executable(linear_nussinov_test
      ../tests/linear_nussinov_test.cpp nussinov.cpp)
  target_include_directories(linear_nussinov_test PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(linear_nussinov_test PRIVATE cxx_std_17)
  add_test(NAME linear_nussinov COMMAND linear_nussinov_test)

  add_executable(ribosum_test ../tests/ribosum_test.cpp ribosum.cpp)
  target_include_directories(ribosum_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
  target_compile_features(ribosum_test PRIVATE cxx_std_17)
  add_test(NAME ribosum_pair_match COMMAND ribosum_test)

  add_test(
      NAME benchmark_scorer
      COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/tests/test_benchmark_scorer.py)

  add_test(
      NAME benchmark_runner_smoke
      COMMAND ${Python3_EXECUTABLE}
        ${PROJECT_SOURCE_DIR}/benchmarks/run.py
        --config ${PROJECT_SOURCE_DIR}/benchmarks/configs/smoke.json
        --binary $<TARGET_FILE:dafs>
        --output-dir ${CMAKE_CURRENT_BINARY_DIR}/benchmark-smoke
        --jobs 2
        --rerun)
  set_tests_properties(benchmark_runner_smoke PROPERTIES TIMEOUT 120)

  add_test(
      NAME reject_negative_max_iter
      COMMAND ${CMAKE_COMMAND}
        -DDAFS=$<TARGET_FILE:dafs>
        -DINPUT=${PROJECT_SOURCE_DIR}/tests/data/tiny.fa
        -P ${PROJECT_SOURCE_DIR}/tests/reject_negative_max_iter.cmake)

  add_test(
      NAME linear_pipeline_smoke
      COMMAND ${CMAKE_COMMAND}
        -DDAFS=$<TARGET_FILE:dafs>
        -DINPUT=${PROJECT_SOURCE_DIR}/tests/data/tiny.fa
        -P ${PROJECT_SOURCE_DIR}/tests/linear_pipeline_smoke.cmake)
  set_tests_properties(linear_pipeline_smoke PROPERTIES TIMEOUT 60)

  if(GLPK_FOUND OR NOT DAFS_HAS_IP_BACKEND)
    add_executable(linear_ipknot_test
        ../tests/linear_ipknot_test.cpp ipknot.cpp ip.cpp nussinov.cpp)
    target_include_directories(linear_ipknot_test PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_features(linear_ipknot_test PRIVATE cxx_std_17)
    if(GLPK_FOUND)
      target_include_directories(linear_ipknot_test PRIVATE ${GLPK_INCLUDE_DIR})
      target_compile_definitions(linear_ipknot_test PRIVATE WITH_GLPK)
      target_link_libraries(linear_ipknot_test PRIVATE ${GLPK_LIBRARY})
    endif()
    add_test(NAME linear_ipknot COMMAND linear_ipknot_test)
  endif()
  add_test(NAME linear_ipknot_pipeline
      COMMAND ${CMAKE_COMMAND}
        -DDAFS=$<TARGET_FILE:dafs>
        -DINPUT=${PROJECT_SOURCE_DIR}/tests/data/tiny.fa
        -P ${PROJECT_SOURCE_DIR}/tests/linear_ipknot_pipeline.cmake)
  set_tests_properties(linear_ipknot_pipeline PROPERTIES TIMEOUT 60)

  if(DAFS_HAS_IP_BACKEND)
    if(GLPK_FOUND)
      add_executable(ipknot_sparse_support_test
          ../tests/ipknot_sparse_support_test.cpp ipknot.cpp ip.cpp nussinov.cpp)
      target_include_directories(ipknot_sparse_support_test PRIVATE
          ${CMAKE_CURRENT_SOURCE_DIR} ${GLPK_INCLUDE_DIR})
      target_compile_features(ipknot_sparse_support_test PRIVATE cxx_std_17)
      target_compile_definitions(ipknot_sparse_support_test PRIVATE WITH_GLPK)
      target_link_libraries(ipknot_sparse_support_test PRIVATE ${GLPK_LIBRARY})
      add_test(NAME ipknot_sparse_support COMMAND ipknot_sparse_support_test)
    endif()
    add_test(
        NAME ipknot_sparse_pipeline
        COMMAND ${CMAKE_COMMAND}
          -DDAFS=$<TARGET_FILE:dafs>
          -DINPUT=${PROJECT_SOURCE_DIR}/tests/data/tiny.fa
          -P ${PROJECT_SOURCE_DIR}/tests/ipknot_sparse_pipeline.cmake)
    set_tests_properties(ipknot_sparse_pipeline PROPERTIES TIMEOUT 60)
    add_test(
        NAME zero_max_iter_uses_ip
        COMMAND ${CMAKE_COMMAND}
          -DDAFS=$<TARGET_FILE:dafs>
          -DINPUT=${PROJECT_SOURCE_DIR}/tests/data/tiny.fa
          -P ${PROJECT_SOURCE_DIR}/tests/zero_max_iter_uses_ip.cmake)
    set_tests_properties(zero_max_iter_uses_ip PROPERTIES TIMEOUT 60)
  endif()
endif()

install(TARGETS dafs RUNTIME DESTINATION bin)
