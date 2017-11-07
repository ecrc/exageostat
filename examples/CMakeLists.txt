set( tests_files zgen_mle_test.c )


foreach(test_src ${tests_files})
    get_filename_component(test_exe ${test_src} NAME_WE)
    add_executable(test_${test_exe} ${test_src})
    target_link_libraries(test_${test_exe} exageostat)
    set_target_properties(test_${test_exe} PROPERTIES OUTPUT_NAME ${test_exe})
endforeach()
