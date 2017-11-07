###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2017 Inria. All rights reserved.
# @copyright (c) 2012-2014 Bordemisc INP, CNRS (LaBRI UMR 5800), Inria, Univ. Bordemisc. All rights reserved.
#
###
#
#  @file CMakeLists.txt
#
#  @project MORSE
#  MORSE is a software package provided by:
#     Inria Bordemisc - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 0.9.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###


# Generate the morse headers for all possible precisions
# ------------------------------------------------------
set(EXAGEOSTAT_HDRS_GENERATED "")
set(ZHDR
#    morse_z.h
#    morse_zc.h
#    runtime_z.h
#    runtime_zc.h
)

precisions_rules_py(EXAGEOSTAT_HDRS_GENERATED "${ZHDR}"
                    PRECISIONS "z;c;d;s;zc;ds")

configure_file("exageostat_config.h.in"
               "exageostat_config.h"
               @ONLY)
configure_file("morse.h.in"
               "morse.h"
               @ONLY)

# Define the list of headers
# --------------------------
set(EXAGEOSTAT_HDRS
    exageostat_config.h
    morse_constants.h
    morse_kernels.h
    morse_simulate.h
    morse_struct.h
    morse_types.h
    morse.h
    runtime.h
    ${EXAGEOSTAT_HDRS_GENERATED}
    )

set(flags_to_add "")
foreach(_prec ${EXAGEOSTAT_PRECISION})
    set(flags_to_add "${flags_to_add} -DPRECISION_${_prec}")
endforeach()
set_source_files_properties(runtime.h PROPERTIES COMPILE_FLAGS "${flags_to_add}")

# Force generation of headers
# ---------------------------
add_custom_target(exageostat_include ALL SOURCES ${EXAGEOSTAT_HDRS})

set(HDR_INSTALL
    morse_constants.h
    morse_fortran.h
    morse_kernels.h
    morse_simulate.h
    morse_struct.h
    morse_types.h
    runtime.h
    )
foreach( hdr_file ${EXAGEOSTAT_HDRS_GENERATED} )
    list(APPEND HDR_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach()
list(APPEND HDR_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/morse.h ${CMAKE_CURRENT_BINARY_DIR}/exageostat_config.h)

# installation
# ------------
install(FILES ${HDR_INSTALL}
        DESTINATION include/exageostat)

###
### END CMakeLists.txt
###
