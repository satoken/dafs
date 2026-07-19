execute_process(
  COMMAND "${DAFS}" -a LinearAlign -s lpc --no-alifold --max-iter=0 "${INPUT}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error)

if(NOT result EQUAL 0)
  message(FATAL_ERROR "--max-iter=0 failed\n${output}${error}")
endif()
if(NOT output MATCHES ">SS_cons[\n]+[.()]+")
  message(FATAL_ERROR "--max-iter=0 produced no consensus structure\n${output}${error}")
endif()
