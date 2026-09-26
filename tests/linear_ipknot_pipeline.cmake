get_filename_component(binary_dir "${DAFS}" DIRECTORY)
set(metrics "${binary_dir}/linear_ipknot_pipeline.metrics.jsonl")

execute_process(
  COMMAND "${DAFS}" -a LinearAlign -s lpc --no-alifold --ipknot
    --max-iter=1 --metrics-jsonl "${metrics}" "${INPUT}"
  RESULT_VARIABLE status
  OUTPUT_VARIABLE prediction
  ERROR_VARIABLE errors)
if(NOT status EQUAL 0)
  message(FATAL_ERROR "Linear IPknot pipeline failed: ${errors}")
endif()
if(NOT prediction MATCHES ">SS_cons[\r\n]+")
  message(FATAL_ERROR "Linear IPknot prediction is missing: ${prediction}")
endif()

file(READ "${metrics}" trace)
if(NOT trace MATCHES "\"structure_decoder\":\"linear-IPknot\"")
  message(FATAL_ERROR "Linear IPknot decoder was not selected")
endif()
if(NOT trace MATCHES "\"sparse_structure_lagrangian\":true")
  message(FATAL_ERROR "Linear IPknot did not use sparse structure storage")
endif()
if(NOT trace MATCHES "\"x_bound\":[0-9]" OR
   NOT trace MATCHES "\"y_bound\":[0-9]")
  message(FATAL_ERROR "Linear IPknot capacity bounds are missing")
endif()
