include_directories("${CMAKE_CURRENT_BINARY_DIR}/../include")

file(GLOB libFiles RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}  "${CMAKE_CURRENT_SOURCE_DIR}/*.cpp") # select all .cpp-files

set(multi_idx_SRCS ${libFiles} )

add_library(multi_idx ${multi_idx_SRCS} )

install(TARGETS multi_idx RUNTIME DESTINATION bin
                    LIBRARY DESTINATION lib
                    ARCHIVE DESTINATION lib)    

set_target_properties(multi_idx PROPERTIES
                      VERSION "1.0.0"   
                      SOVERSION "1"    
)
