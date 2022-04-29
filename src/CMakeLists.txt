set(FQFEEDER_LIB_SRCS ${FQFEEDER_LIB_DIR}/FastxParser.cpp)

add_library(fqfeeder ${FQFEEDER_LIB_SRCS})
target_include_directories(fqfeeder PUBLIC ${FQFEEDER_INCLUDE_DIR})

add_executable(generate_count_matrix_ADTs generate_count_matrix_ADTs.cpp) 
target_include_directories(generate_count_matrix_ADTs PUBLIC ${FQFEEDER_INCLUDE_DIR})
target_link_libraries(generate_count_matrix_ADTs pthread fqfeeder z)
