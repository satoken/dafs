execute_process(
  COMMAND "${DAFS}" --max-iter=-1 "${INPUT}"
  RESULT_VARIABLE result
  OUTPUT_VARIABLE output
  ERROR_VARIABLE error)

if(result EQUAL 0)
  message(FATAL_ERROR "negative --max-iter unexpectedly succeeded\n${output}${error}")
endif()
if(NOT error MATCHES "--max-iter must be non-negative")
  message(FATAL_ERROR "missing max-iter diagnostic\n${output}${error}")
endif()
