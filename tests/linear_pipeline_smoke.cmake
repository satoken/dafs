execute_process(
  COMMAND "${DAFS}" -a LinearAlign -s lpc --no-alifold
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
