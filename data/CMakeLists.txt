SET(base_url "http://xor.di.unipi.it/~rossano")

# A samll default data set for testing
SET(local_test_cases test.hash)
## Uncomment the following line to use the datasets from the SIGIR 2016 paper
#SET(local_test_cases Clueweb09-Full.SimHash Clueweb09-Full.OddSketch lsh_sift_64.hash mlh_sift_64.hash)
SET(test_cases ${local_test_cases} PARENT_SCOPE)

FOREACH(test_case ${local_test_cases})
    SET(abs_test_case ${CMAKE_HOME_DIRECTORY}/data/${test_case})
    IF(NOT EXISTS ${abs_test_case})
        SET(tc_url ${base_url}/${test_case})
        MESSAGE(${tc_url})
        FILE(DOWNLOAD ${tc_url} ${abs_test_case} SHOW_PROGRESS)
    endif(NOT EXISTS ${abs_test_case})
ENDFOREACH(test_case)
