execute_process(
  COMMAND "${DAFS}" -a LinearAlign -s lpc
          --align-beam=7 --align-dd-beam=9
          --linfold-beam=11 --fold-dd-beam=13 --fold-final-beam=15
          --ribosum-weight=0.1 --max-iter=2 "${INPUT}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error)

if(NOT result EQUAL 0)
  message(FATAL_ERROR "linear pipeline failed\n${output}${error}")
endif()
if(NOT output MATCHES ">SS_cons[\n]+[.()]+")
  message(FATAL_ERROR "linear pipeline produced no consensus structure\n${output}${error}")
endif()
if(error MATCHES "LinearPartition.*failed|LinearAlign failed")
  message(FATAL_ERROR "linear model failure was reported\n${output}${error}")
endif()

execute_process(
  COMMAND "${DAFS}" -a LinearAlign -s lpc --no-alifold
          --align-beam=7 --align-dd-beam=9
          --linfold-beam=11 --fold-dd-beam=13 --fold-final-beam=15
          --ribosum-weight=0.1 --max-iter=2 "${INPUT}"
  RESULT_VARIABLE no_alifold_result
  OUTPUT_VARIABLE no_alifold_output
  ERROR_VARIABLE no_alifold_error)
if(NOT no_alifold_result EQUAL 0)
  message(FATAL_ERROR
          "linear --no-alifold pipeline failed\n${no_alifold_output}${no_alifold_error}")
endif()
if(NOT output STREQUAL no_alifold_output)
  message(FATAL_ERROR
          "linear profile surrogate was not disabled consistently")
endif()

# The public defaults are RIBOSUM 0.075 and RNAalifold disabled.  Check the
# nonlinear path because the linear folding engine disables RNAalifold on its
# own to preserve linear complexity.
execute_process(
  COMMAND "${DAFS}" -a CONTRAlign -s CONTRAfold --dynamic-cbp --max-iter=2
          "${INPUT}"
  RESULT_VARIABLE default_result
  OUTPUT_VARIABLE default_output
  ERROR_VARIABLE default_error)
if(NOT default_result EQUAL 0)
  message(FATAL_ERROR
          "default-option pipeline failed\n${default_output}${default_error}")
endif()

execute_process(
  COMMAND "${DAFS}" -a CONTRAlign -s CONTRAfold --no-alifold --dynamic-cbp
          --ribosum-weight=0.075 --max-iter=2 "${INPUT}"
  RESULT_VARIABLE explicit_default_result
  OUTPUT_VARIABLE explicit_default_output
  ERROR_VARIABLE explicit_default_error)
if(NOT explicit_default_result EQUAL 0)
  message(FATAL_ERROR
          "explicit-default pipeline failed\n${explicit_default_output}${explicit_default_error}")
endif()
if(NOT default_output STREQUAL explicit_default_output)
  message(FATAL_ERROR
          "implicit defaults differ from --ribosum-weight=0.075 --no-alifold")
endif()

execute_process(
  COMMAND "${DAFS}" --alifold --no-alifold "${INPUT}"
  RESULT_VARIABLE conflicting_alifold_result
  OUTPUT_QUIET
  ERROR_QUIET)
if(conflicting_alifold_result EQUAL 0)
  message(FATAL_ERROR "--alifold and --no-alifold were accepted together")
endif()
