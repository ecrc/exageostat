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
#  @project CHAMELEON
#  CHAMELEON is a software package provided by:
#     Inria Bordemisc - Sud-Ouest,
#     Univ. of Tennessee,
#     King Abdullah Univesity of Science and Technology
#     Univ. of California Berkeley,
#     Univ. of Colorado Denver.
#
#  @version 1.2.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###


# Generate the chameleon headers for all possible precisions
# ------------------------------------------------------
set(EXAGEOSTAT_HDRS_GENERATED "")
set(ZHDR)

precisions_rules_py(EXAGEOSTAT_HDRS_GENERATED "${ZHDR}"
        PRECISIONS "z;c;d;s;zc;ds")

configure_file("exageostat_config.h.in"
        "exageostat_config.h"
        @ONLY)
configure_file("chameleon/chameleon_config.h.in"
        "chameleon/chameleon_config.h"
        @ONLY)


# Define the list of headers
# --------------------------
set(EXAGEOSTAT_HDRS
        exageostat_config.h
        constants.h
        simulate.h
        struct.h
        tasks.h
        types.h
        chameleon.h
        runtime.h
        ${EXAGEOSTAT_HDRS_GENERATED}
        )

set(flags_to_add "")
foreach (_prec ${EXAGEOSTAT_PRECISION})
    set(flags_to_add "${flags_to_add} -DPRECISION_${_prec}")
endforeach ()
set_source_files_properties(runtime.h PROPERTIES COMPILE_FLAGS "${flags_to_add}")

# Force generation of headers
# ---------------------------
add_custom_target(exageostat_include ALL SOURCES ${EXAGEOSTAT_HDRS})

set(HDR_INSTALL
        constants.h
        fortran.h
        struct.h
        types.h
        runtime.h
        )
foreach (hdr_file ${EXAGEOSTAT_HDRS_GENERATED})
    list(APPEND HDR_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/${hdr_file})
endforeach ()
list(APPEND HDR_INSTALL ${CMAKE_CURRENT_BINARY_DIR}/chameleon.h ${CMAKE_CURRENT_BINARY_DIR}/exageostat_config.h)

# installation
# ------------
install(FILES ${HDR_INSTALL}
        DESTINATION include/exageostat)

###
### END CMakeLists.txt
###
