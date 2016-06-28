cmake_minimum_required(VERSION 2.8)
cmake_policy(SET CMP0015 NEW)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(AppendCompilerFlags)
include(ExternalProject)
include(MultiIndexHelper)

enable_testing()

project(multi_idx CXX C)

append_cxx_compiler_flags("-msse4.2 -std=c++14 -Wall -DNDEBUG" "GCC" CMAKE_CXX_FLAGS)
append_cxx_compiler_flags("-O3 -ffast-math -funroll-loops" "GCC" CMAKE_CXX_FLAGS)



LINK_DIRECTORIES(${CMAKE_BINARY_DIR}/external/sdsl-lite/lib)
INCLUDE_DIRECTORIES(${CMAKE_HOME_DIRECTORY}/include
                    ${CMAKE_BINARY_DIR}/external/sdsl-lite/include
                    ${CMAKE_BINARY_DIR}/external/sdsl-lite/external/libdivsufsort-2.0.1/include
                    )

find_package(Git)
if(GIT_FOUND)
    message("git found: ${GIT_EXECUTABLE}")
else()
    message(WARNING "git not found. Cloning of submodules will not work.")
endif()

ADD_SUBDIRECTORY(external/sdsl-lite)
ADD_SUBDIRECTORY(lib)

ADD_EXECUTABLE(gen_bench_data src/gen_bench_data.cpp)
TARGET_LINK_LIBRARIES(gen_bench_data sdsl)

ADD_EXECUTABLE(bench_scan src/bench_scan.cpp)
TARGET_LINK_LIBRARIES(bench_scan sdsl)

ADD_EXECUTABLE(bench_scan2 src/bench_scan2.cpp)
TARGET_LINK_LIBRARIES(bench_scan2 sdsl multi_idx)

ADD_EXECUTABLE(bench_scan3 src/bench_scan3.cpp)
TARGET_LINK_LIBRARIES(bench_scan3 sdsl multi_idx)

#ADD_EXECUTABLE(cluster_statistics src/cluster_statistics.cpp)
#TARGET_LINK_LIBRARIES(cluster_statistics sdsl multi_idx)

ADD_SUBDIRECTORY(data)

#  Generate target for the construction of key databases and queries
FOREACH(test_case ${test_cases})
    SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
        ADD_CUSTOM_COMMAND(OUTPUT ${abs_test_case}.data
                                  ${abs_test_case}.query
                                  ${abs_test_case}.100.query
                           COMMAND $<TARGET_FILE:gen_bench_data> ${abs_test_case} 10000
                           COMMAND cat ${abs_test_case}.existing.query ${abs_test_case}.real.query > ${abs_test_case}.query
                           COMMAND head -c 400 ${abs_test_case}.existing.query > ${abs_test_case}.100.query
                           COMMAND head -c 400 ${abs_test_case}.real.query >> ${abs_test_case}.100.query
                           DEPENDS gen_bench_data
                           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                           COMMENT "Generate query files for ${test_case}."
                           VERBATIM)
        ADD_CUSTOM_TARGET(${test_case}-query-files
                          DEPENDS ${abs_test_case}.data
                                  ${abs_test_case}.query)
        MESSAGE("Add target ${test_case}-query-files")
ENDFOREACH(test_case)

SET(perm_header_cmake "${CMAKE_HOME_DIRECTORY}/include/multi_idx/perm.hpp.cmake")
SET(perm_header "${CMAKE_HOME_DIRECTORY}/include/multi_idx/perm.hpp")
FILE(STRINGS ${perm_header_cmake} header NEWLINE_CONSUME)
FILE(REMOVE ${perm_header})
FILE(APPEND ${perm_header} ${header})

ADD_EXECUTABLE(gen_hash_file src/gen_hash_file.cpp)
TARGET_LINK_LIBRARIES(gen_hash_file sdsl)

ADD_EXECUTABLE(sim_hash src/sim_hash.cpp)

## Input info

ADD_EXECUTABLE(ham_distribution src/ham_distribution)
TARGET_LINK_LIBRARIES(ham_distribution sdsl)

SET(input_infos "")

FOREACH(test_case ${test_cases})
    SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case}.data)
    SET(input_info ${abs_test_case}.distr)
    LIST(APPEND input_infos ${input_info})
    ADD_CUSTOM_COMMAND(OUTPUT ${input_info}
                       COMMAND $<TARGET_FILE:ham_distribution> ${abs_test_case} > ${input_info}
                       DEPENDS ham_distribution
                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                       COMMENT "Calculating hamming weight distribution for ${test_case}" VERBATIM)
ENDFOREACH(test_case) 

ADD_CUSTOM_TARGET(input_info DEPENDS ${input_infos})


## Experiment 0: Compare mi and mi_red
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/results)
SET(exp0_result_file "${CMAKE_BINARY_DIR}/results/exp0.result.txt")
SET(exp0_results "")

FILE(STRINGS ${CMAKE_HOME_DIRECTORY}/exp0.config exp0_lines REGEX "^[^#].*")
FOREACH(line ${exp0_lines})
    MESSAGE("line = ${line}")
    LIST(GET line 0 index_name)
    LIST(GET line 1 index_type)
    LIST(GET line 2 errors)
    STRING(REPLACE "," ";" error_list ${errors})
    MESSAGE("Exp 0 ${index_name}; ${index_type}; ${error_list}")
    FOREACH(t_k ${error_list})
        STRING(REGEX REPLACE "([,<])t_k([,>])" "\\1${t_k}\\2" type ${index_type}) 
        GEN_PERM_FILE(${type} blocks)
        SET(exec ${index_name}_index_${t_k})
        IF(NOT TARGET ${exec})
            ADD_EXECUTABLE(${exec} src/index.cpp)
            TARGET_LINK_LIBRARIES(${exec} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
            SET_PROPERTY(TARGET ${exec} PROPERTY COMPILE_DEFINITIONS 
                         BLOCKS=${blocks}
                         K=${t_k} 
                         INDEX_TYPE=${index_type} 
                         INDEX_NAME="${index_name}")

            FOREACH(test_case ${test_cases})
                SET(exp0_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.query.exp0.result.txt)
                LIST(APPEND exp0_results ${exp0_result})
                SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
                IF(NOT TARGET ${exp0_result})
                    ADD_CUSTOM_COMMAND(OUTPUT ${exp0_result}
#                                   COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
#                                   COMMAND sh -c "sync && purge"                        
                                       COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.query 0 > ${exp0_result}
#                                   COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
                                       COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.query 1 >> ${exp0_result}
                                       COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.100.query 0 1
                                       DEPENDS ${exec} ${test_case}-query-files
                                       WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                                       COMMENT "Execute exp0 for t_k=${t_k}\nCreating ${exp0_result}.\n"
                                       VERBATIM)
                ENDIF()
            ENDFOREACH()
        ENDIF()
    ENDFOREACH()
ENDFOREACH()

ADD_CUSTOM_TARGET(exp0
                  COMMAND cat ${exp0_results} > ${exp0_result_file}
                  COMMAND R CMD BATCH ../scripts/exp0.R
                  DEPENDS ${exp0_results}
                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
                  COMMENT "\nResults can be found in ${exp0_result_file}.\nVisualization in file exp0.pdf\n")

### Experiment 1: Measure the efficiency of the searching phase
#SET(exp1_result_file "${CMAKE_BINARY_DIR}/results/exp1.result.txt")
#SET(exp1_results "")
#
#FILE(STRINGS ${CMAKE_HOME_DIRECTORY}/exp1.config exp1_lines REGEX "^[^#].*")
#FOREACH(line ${exp1_lines})
#    LIST(GET line 0 index_name)
#    LIST(GET line 1 index_type)
#    LIST(GET line 2 blocks)
#    LIST(GET line 3 errors)
#    MESSAGE("Exp 1 ${index_name}; ${index_type}; ${blocks}; ${errors}")
#
#    GEN_PERM_FILE(${blocks} ${errors} ${index_type})
#
#    SET(exec ${index_name}_index_${blocks}_${errors})
#    IF(NOT TARGET ${exec})
#        ADD_EXECUTABLE(${exec} src/index.cpp)
#        TARGET_LINK_LIBRARIES(${index_name}_index_${blocks}_${errors} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
#        SET_PROPERTY(TARGET ${index_name}_index_${blocks}_${errors} PROPERTY COMPILE_DEFINITIONS 
#                     N=${blocks}
#                     K=${errors} 
#                     INDEX_TYPE=${index_type} 
#                     INDEX_NAME="${index_name}")
#    ENDIF(NOT TARGET ${exec})
#
#
#    FOREACH(test_case ${test_cases})
#        SET(queries query)
#        FOREACH(query_suf ${queries})
#            SET(exp1_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.${query_suf}.exp1.result.txt)
#            LIST(APPEND exp1_results ${exp1_result})
#            SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
#            IF(NOT TARGET ${exp1_result})
#                ADD_CUSTOM_COMMAND(OUTPUT ${exp1_result}
#                                   COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
#                                   COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.${query_suf} 0 > ${exp1_result}
#                                   COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
#                                   COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.${query_suf} 1 >> ${exp1_result}
#                                   DEPENDS ${exec} ${test_case}-query-files
#                                   WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                                   COMMENT "Execute exp1 for b=${blocks} and e=${errors}\nCreating ${exp1_result}.\n \
#                                            $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.${query_suf} > ${exp1_result}"
#                                   VERBATIM)
#
#            ENDIF(NOT TARGET ${exp1_result})
#        ENDFOREACH(query_suf)
#    ENDFOREACH(test_case)
#ENDFOREACH(line)
#
#ADD_CUSTOM_TARGET(exp1
#                  COMMAND cat ${exp1_results} > ${exp1_result_file}
##                  COMMAND R --vanilla < ../scripts/table-1.R
#                  DEPENDS ${exp1_results}
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                  COMMENT "Summarize results in ${exp1_result_file}.")
#
## Experiment 2
#
#SET(exp2_result_file "${CMAKE_BINARY_DIR}/results/exp2.result.txt")
#SET(exp2_results "")
#
#
#FILE(STRINGS ${CMAKE_HOME_DIRECTORY}/exp2.config exp2_lines REGEX "^[^#].*")
#FOREACH(line ${exp2_lines})
#    LIST(GET line 0 index_name)
#    LIST(GET line 1 index_type)
#    LIST(GET line 2 blocks)
#    LIST(GET line 3 errors)
#    LIST(GET line 4 query)
#    SET(query_suf ${query}.query)
#    MESSAGE("Exp 2 ${index_name}; ${index_type}; ${blocks}; ${errors}")
#
#    GEN_PERM_FILE(${blocks} ${errors} ${index_type})
#
#    SET(exec ${index_name}_index_${blocks}_${errors})
#    IF(NOT TARGET ${exec})
#        MESSAGE("Add target ${exec}")
#        ADD_EXECUTABLE(${exec} src/index.cpp)
#        TARGET_LINK_LIBRARIES(${index_name}_index_${blocks}_${errors} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
#        SET_PROPERTY(TARGET ${index_name}_index_${blocks}_${errors} PROPERTY COMPILE_DEFINITIONS 
#                     N=${blocks}
#                     K=${errors} 
#                     INDEX_TYPE=${index_type} 
#                     INDEX_NAME="${index_name}")
#    ENDIF(NOT TARGET ${exec})
#
#
#    FOREACH(test_case ${test_cases})
#        SET(exp2_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.${query_suf}.exp2.result.txt)
#        LIST(APPEND exp2_results ${exp2_result})
#        SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
#        IF(NOT TARGET ${exp2_result})
#            MESSAGE("Add command ${exp2_result}")
#            ADD_CUSTOM_COMMAND(OUTPUT ${exp2_result}
#                               COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
#                               COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.${query_suf} 0  > ${exp2_result}
#                               COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
#                               COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.${query_suf} 1 >> ${exp2_result}
#                               DEPENDS ${exec} ${test_case}-query-files
#                               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                               COMMENT "Execute exp2 for b=${blocks} and e=${errors}\nCreating ${exp2_result}.\n \
#                                        $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.${query_suf} > ${exp2_result}"
#                               VERBATIM)
#        ENDIF(NOT TARGET ${exp2_result})
#    ENDFOREACH(test_case)
#ENDFOREACH(line)
#
#ADD_CUSTOM_TARGET(exp2
#                  COMMAND cat ${exp2_results} > ${exp2_result_file}
#                  DEPENDS ${exp2_results}
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                  COMMENT "Summarize results in ${exp2_result_file}.")
#
### Experiment 3: Evaluate the effect of splitting on mi_bv_red
#
#SET(exp3_result_file "${CMAKE_BINARY_DIR}/results/exp3.result.txt")
#SET(exp3_results "")
#
#FILE(STRINGS ${CMAKE_HOME_DIRECTORY}/exp3.config exp3_lines REGEX "^[^#].*")
#FOREACH(line ${exp3_lines})
#    LIST(GET line 0 index_name)
#    LIST(GET line 1 index_type)
#    LIST(GET line 2 blocks)
#    LIST(GET line 3 errors)
#    MESSAGE("Exp 3 ${index_name}; ${index_type}; ${blocks}; ${errors}")
#
#    GEN_PERM_FILE(${blocks} ${errors} ${index_type})
#
#    SET(exec ${index_name}_index_${blocks}_${errors})
#    IF(NOT TARGET ${exec})
#        ADD_EXECUTABLE(${exec} src/index.cpp)
#        TARGET_LINK_LIBRARIES(${index_name}_index_${blocks}_${errors} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
#        SET_PROPERTY(TARGET ${index_name}_index_${blocks}_${errors} PROPERTY COMPILE_DEFINITIONS 
#                     N=${blocks}
#                     K=${errors} 
#                     INDEX_TYPE=${index_type} 
#                     INDEX_NAME="${index_name}")
#    ENDIF(NOT TARGET ${exec})
#
#
#    SET(test_case "lsh_sift_64.hash")
#    SET(exp3_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.exp3.result.txt)
#    LIST(APPEND exp3_results ${exp3_result})
#    SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
#
#    SET(input_lines "")
#    FOREACH(query_suf ${queries})
#        SET(input_lines "${abs_test_case}.${query_suf}:${input_lines}")
#    ENDFOREACH(query_suf)
#
#    IF(NOT TARGET ${exp3_result})
#        ADD_CUSTOM_COMMAND(OUTPUT ${exp3_result}
#                           COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
#                           COMMAND echo ${input_lines} | $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${exp3_result}
#                           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                           COMMENT "Execute exp3 for b=${blocks} and e=${errors}\nCreating ${exp3_result}.\n \
#                                    echo ${input_lines} | $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${exp3_result}"
#                           VERBATIM)
#
#    ENDIF(NOT TARGET ${exp3_result})
#ENDFOREACH(line)
#
#ADD_CUSTOM_TARGET(exp3
#                  COMMAND cat ${exp3_results} > ${exp3_result_file}
##                  COMMAND R --vanilla < ../scripts/table-3.R
#                  DEPENDS ${exp3_results}
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                  COMMENT "Summarize results in ${exp3_result_file}.")
#
#
#
### Experiment 5: Evaluate the best threshold for triangle cluster
#SET(exp5_result_file "${CMAKE_BINARY_DIR}/results/exp5.result.txt")
#SET(exp5_results "")
#SET(expRS_result_file "${CMAKE_BINARY_DIR}/results/expRS.result.txt")
#SET(expRS_results "")
#SET(queries query)
#
#FILE(STRINGS ${CMAKE_HOME_DIRECTORY}/exp5.config exp5_lines REGEX "^[^#].*")
#FOREACH(line ${exp5_lines})
#    LIST(GET line 0 index_name)
#    LIST(GET line 1 index_type)
#    LIST(GET line 2 blocks)
#    LIST(GET line 3 errors)
#    MESSAGE("Exp 5 ${index_name}; ${index_type}; ${blocks}; ${errors}")
#
#    GEN_PERM_FILE(${blocks} ${errors} ${index_type})
#
#    SET(exec ${index_name}_index_${blocks}_${errors})
#    SET(exec_rs ${index_name}_index_rs_${blocks}_${errors})
#    IF(NOT TARGET ${exec})
#        ADD_EXECUTABLE(${exec} src/index.cpp)
#        TARGET_LINK_LIBRARIES(${index_name}_index_${blocks}_${errors} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
#        SET_PROPERTY(TARGET ${index_name}_index_${blocks}_${errors} PROPERTY COMPILE_DEFINITIONS 
#                     N=${blocks}
#                     K=${errors} 
#                     INDEX_TYPE=${index_type} 
#                     INDEX_NAME="${index_name}")
#    ENDIF(NOT TARGET ${exec})
#    ADD_EXECUTABLE(${exec_rs} src/index.cpp)
#    TARGET_LINK_LIBRARIES(${index_name}_index_rs_${blocks}_${errors} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
#    SET_PROPERTY(TARGET ${index_name}_index_rs_${blocks}_${errors} PROPERTY COMPILE_DEFINITIONS 
#                 N=${blocks}
#                 K=${errors} 
#                 STATS=1
#                 INDEX_TYPE=${index_type} 
#                 INDEX_NAME="${index_name}")
#
#
#
#    FOREACH(test_case ${test_cases})
#        SET(exp5_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.exp5.result.txt)
#        SET(expRS_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.expRS.result.txt)
#        LIST(APPEND exp5_results ${exp5_result})
#        LIST(APPEND expRS_results ${expRS_result})
#        SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
#
#        SET(input_lines "")
#        FOREACH(query_suf ${queries})
#            SET(input_lines "${abs_test_case}.${query_suf} 0 0 0 1:${input_lines}")
#        ENDFOREACH(query_suf)
#
#        IF(NOT TARGET ${exp5_result})
#            ADD_CUSTOM_COMMAND(OUTPUT ${exp5_result}
#                               COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
#                               COMMAND echo ${input_lines} | $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${exp5_result}
#                               DEPENDS ${exec} ${test_case}-query-files
#                               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                               COMMENT "Execute exp5 for b=${blocks} and e=${errors}\nCreating ${exp5_result}.\n \
#                                        echo ${input_lines} | $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${exp5_result}"
#                               VERBATIM)
#        ENDIF(NOT TARGET ${exp5_result})
#
#        IF(NOT TARGET ${expRS_result})
#            ADD_CUSTOM_COMMAND(OUTPUT ${expRS_result}
#                               COMMAND echo ${input_lines} | $<TARGET_FILE:${exec_rs}> ${abs_test_case}.data > ${expRS_result}
#                               DEPENDS ${exec_rs} ${test_case}-query-files
#                               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                               COMMENT "Execute expRS for b=${blocks} and e=${errors}\nCreating ${expRS_result}.\n \
#                                        echo ${input_lines} | $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${expRS_result}"
#                               VERBATIM)
#        ENDIF(NOT TARGET ${expRS_result})
#
#    ENDFOREACH(test_case)
#ENDFOREACH(line)
#
#ADD_CUSTOM_TARGET(exp5
#                  COMMAND cat ${exp5_results} > ${exp5_result_file}
#                  DEPENDS ${exp5_results}
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                  COMMENT "Summarize results in ${exp5_result_file}.")
#
#ADD_CUSTOM_TARGET(expRS
#                  COMMAND cat ${expRS_results} > ${expRS_result_file}
#                  DEPENDS ${expRS_results}
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                  COMMENT "Summarize results in ${expRS_result_file}.")
#
#
#
### Experiment 6: Final Comparison
#SET(exp6_result_file "${CMAKE_BINARY_DIR}/results/exp6.result.txt")
#SET(exp6_results "")
#SET(queries query)
#
#FILE(STRINGS ${CMAKE_HOME_DIRECTORY}/exp6.config exp6_lines REGEX "^[^#].*")
#FOREACH(line ${exp6_lines})
#    LIST(GET line 0 index_name)
#    LIST(GET line 1 index_type)
#    LIST(GET line 2 blocks)
#    LIST(GET line 3 errors)
#    MESSAGE("Exp 6 ${index_name}; ${index_type}; ${blocks}; ${errors}")
#
#    GEN_PERM_FILE(${blocks} ${errors} ${index_type})
#
#    SET(exec ${index_name}_index_${blocks}_${errors})
#    IF(NOT TARGET ${exec})
#        ADD_EXECUTABLE(${exec} src/index.cpp)
#        TARGET_LINK_LIBRARIES(${index_name}_index_${blocks}_${errors} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
#        SET_PROPERTY(TARGET ${index_name}_index_${blocks}_${errors} PROPERTY COMPILE_DEFINITIONS 
#                     N=${blocks}
#                     K=${errors} 
#                     INDEX_TYPE=${index_type} 
#                     INDEX_NAME="${index_name}")
#    ENDIF(NOT TARGET ${exec})
#
#
#    FOREACH(test_case ${test_cases})
#        SET(exp6_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.exp6.result.txt)
#        LIST(APPEND exp6_results ${exp6_result})
#        SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
#
#        SET(input_lines "")
#        FOREACH(query_suf ${queries})
#            SET(input_lines "${abs_test_case}.${query_suf} 0 0 0 1:${input_lines}")
#        ENDFOREACH(query_suf)
#
#        IF(NOT TARGET ${exp6_result})
#            ADD_CUSTOM_COMMAND(OUTPUT ${exp6_result}
#                               COMMAND sh -c "echo 1 >/proc/sys/vm/drop_caches"
#COMMAND echo ${input_lines} | $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${exp6_result}
##                               COMMAND echo "" | $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${exp6_result}
#                               DEPENDS ${exec} ${test_case}-query-files
#                               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                               COMMENT "Execute exp6 for b=${blocks} and e=${errors}\nCreating ${exp6_result}.\n \
#                                        echo ${input_lines} | $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${exp6_result}"
#                               VERBATIM)
#        ENDIF(NOT TARGET ${exp6_result})
#    ENDFOREACH(test_case)
#ENDFOREACH(line)
#
#ADD_CUSTOM_TARGET(exp6
#                  COMMAND cat ${exp6_results} > ${exp6_result_file}
#                  DEPENDS ${exp6_results}
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                  COMMENT "Summarize results in ${exp6_result_file}.")
#
#
#
### Experiment S: Statistics
#SET(expS_result_file "${CMAKE_BINARY_DIR}/results/expS.result.txt")
#SET(expS_results "")
#SET(queries query)
#
#FILE(STRINGS ${CMAKE_HOME_DIRECTORY}/expS.config expS_lines REGEX "^[^#].*")
#FOREACH(line ${expS_lines})
#    LIST(GET line 0 index_name)
#    LIST(GET line 1 index_type)
#    LIST(GET line 2 blocks)
#    LIST(GET line 3 errors)
#    MESSAGE("Exp S ${index_name}; ${index_type}; ${blocks}; ${errors}")
#    GEN_PERM_FILE(${blocks} ${errors} ${index_type})
#
#    SET(exec ${index_name}_stat_${blocks}_${errors})
#    IF(NOT TARGET ${exec})
#        ADD_EXECUTABLE(${exec} src/cluster_statistics.cpp)
#        TARGET_LINK_LIBRARIES(${index_name}_stat_${blocks}_${errors} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
#        SET_PROPERTY(TARGET ${index_name}_stat_${blocks}_${errors} PROPERTY COMPILE_DEFINITIONS 
#                     N=${blocks}
#                     K=${errors} 
#                     INDEX_TYPE=${index_type} 
#                     INDEX_NAME="${index_name}")
#    ENDIF(NOT TARGET ${exec})
#
#
#    FOREACH(test_case ${test_cases})
#        SET(expS_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.expS.result.txt)
#        LIST(APPEND expS_results ${expS_result})
#        SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
#
#        IF(NOT TARGET ${expS_result})
#            ADD_CUSTOM_COMMAND(OUTPUT ${expS_result}
#                               COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${expS_result}
#                               DEPENDS ${exec} ${abs_test_case}.data
#                               WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                               COMMENT "Execute expS for b=${blocks} and e=${errors}\nCreating ${expS_result}.\n \
#                                        $<TARGET_FILE:${exec}> ${abs_test_case}.data > ${expS_result}"
#                               VERBATIM)
#        ENDIF(NOT TARGET ${expS_result})
#    ENDFOREACH(test_case)
#ENDFOREACH(line)
#
#
#ADD_CUSTOM_TARGET(expS
#                  COMMAND cat ${expS_results} > ${expS_result_file}
#                  DEPENDS ${expS_results}
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                  )
#
#
### Experiment linear scan
#SET(linear_scan_result_file "${CMAKE_BINARY_DIR}/results/linear_scan.result.txt")
#SET(linear_scan_results "")
#SET(queries query)
#
#SET(linear_scan_ks 2;3;4;5;6;7;8;9)
#FOREACH(k ${linear_scan_ks})
#    SET(index_name "linear_scan")
#    SET(index_type "linear_scan<t_k>")
#    SET(errors ${k})
#    MESSAGE("Experiment linear scan ${index_name}; ${index_type}; ${errors}")
#
#    SET(exec ${index_name}_${errors})
#    ADD_EXECUTABLE(${exec} src/index.cpp)
#    TARGET_LINK_LIBRARIES(${exec} sdsl -ggdb divsufsort divsufsort64 multi_idx pthread)
#        SET_PROPERTY(TARGET ${exec} PROPERTY COMPILE_DEFINITIONS 
#                     N=0
#                     K=${errors} 
#                     INDEX_TYPE=${index_type} 
#                     INDEX_NAME="${index_name}")
#
#    FOREACH(test_case ${test_cases})
##GET_FILENAME_COMPONENT(test_case ${test_case} NAME_WT)
#        SET(linear_scan_result ${CMAKE_BINARY_DIR}/results/${exec}.${test_case}.result.txt)
#        LIST(APPEND linear_scan_results ${linear_scan_result})
#        SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
#MESSAGE("ABS_TEST_CASE=${abs_test_case}")
#        ADD_CUSTOM_COMMAND(OUTPUT ${linear_scan_result}
#                           COMMAND head -c 80 ${abs_test_case}.query > ${abs_test_case}.query.10
#                           COMMAND $<TARGET_FILE:${exec}> ${abs_test_case}.data ${abs_test_case}.query.10 > ${linear_scan_result}
#                           DEPENDS ${exec} ${abs_test_case}.data ${abs_test_case}.query
#                           WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                           COMMENT "Execute linear_scan for e=${errors}\nCreating ${linear_scan_result}.\n"
#                           VERBATIM)
#    ENDFOREACH(test_case)
#ENDFOREACH(k)
#
#
#ADD_CUSTOM_TARGET(exp_linear_scan
#                  COMMAND cat ${linear_scan_results} > ${linear_scan_result_file}
#                  DEPENDS ${linear_scan_results}
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
#                  )
#
