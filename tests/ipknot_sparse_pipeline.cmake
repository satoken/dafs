set(common_args
  -a LinearAlign -s lpc --no-alifold --ipknot
  -g 99,99 -G 99,99 -w 20 --max-iter=3)

execute_process(
  COMMAND "${DAFS}" ${common_args} "${INPUT}"
  RESULT_VARIABLE sparse_status
  OUTPUT_VARIABLE sparse_output
  ERROR_VARIABLE sparse_error)
if(NOT sparse_status EQUAL 0)
  message(FATAL_ERROR "Sparse IPknot decoding failed: ${sparse_error}")
endif()

execute_process(
  COMMAND "${DAFS}" ${common_args} --dense-lagrangian "${INPUT}"
  RESULT_VARIABLE dense_status
  OUTPUT_VARIABLE dense_output
  ERROR_VARIABLE dense_error)
if(NOT dense_status EQUAL 0)
  message(FATAL_ERROR "Dense IPknot decoding failed: ${dense_error}")
endif()

if(NOT sparse_output STREQUAL dense_output)
  message(FATAL_ERROR
    "Sparse and dense IPknot predictions differ:\n${sparse_output}\n${dense_output}")
endif()
if(NOT sparse_output MATCHES ">SS_cons[\r\n]+[(][(][(][.][.][.][)][)][)]")
  message(FATAL_ERROR "IPknot test did not select the expected base pairs: ${sparse_output}")
endif()

set(linear_metrics "${DAFS}.linear-ipknot.metrics.jsonl")
file(REMOVE "${linear_metrics}")
execute_process(
  COMMAND "${DAFS}" ${common_args} --metrics-jsonl "${linear_metrics}" "${INPUT}"
  RESULT_VARIABLE linear_status
  OUTPUT_QUIET
  ERROR_VARIABLE linear_error)
if(NOT linear_status EQUAL 0)
  message(FATAL_ERROR "Linear IPknot decoding failed: ${linear_error}")
endif()
file(READ "${linear_metrics}" linear_metric_text)
if(NOT linear_metric_text MATCHES "\"structure_decoder\":\"linear-IPknot\"")
  message(FATAL_ERROR "Linear IPknot decoder was not selected: ${linear_metric_text}")
endif()

execute_process(
  COMMAND "${DAFS}" -a LinearAlign -s lpc --no-alifold
    --ipknot -g 99,99 -G 99,99 -w 20 --max-iter=0
    --metrics-jsonl "${DAFS}.exact-ip.metrics.jsonl" "${INPUT}"
  RESULT_VARIABLE ip_status
  OUTPUT_VARIABLE ip_output
  ERROR_VARIABLE ip_error)
if(NOT ip_status EQUAL 0 OR
   NOT ip_output MATCHES ">SS_cons[\r\n]+[(][(][(][.][.][.][)][)][)]")
  message(FATAL_ERROR "Sparse exact IP path failed: ${ip_output}${ip_error}")
endif()
file(READ "${DAFS}.exact-ip.metrics.jsonl" exact_metric_text)
if(NOT exact_metric_text MATCHES "\"structure_decoder\":\"exact-coupled-IP\"")
  message(FATAL_ERROR "Exact coupled IP path was not reported: ${exact_metric_text}")
endif()
