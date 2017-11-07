configure_file("${CMAKE_CURRENT_SOURCE_DIR}/config.in"
    "${CMAKE_CURRENT_BINARY_DIR}/config")
add_custom_target(docs COMMAND
    "${DOXYGEN_EXECUTABLE}" "${CMAKE_CURRENT_BINARY_DIR}/config")
