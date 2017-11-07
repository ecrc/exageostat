###
#
# @copyright (c) 2009-2014 The University of Tennessee and The University
#                          of Tennessee Research Foundation.
#                          All rights reserved.
# @copyright (c) 2012-2016 Inria. All rights reserved.
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
#  @version 0.9.0
#  @author Cedric Castagnede
#  @author Emmanuel Agullo
#  @author Mathieu Faverge
#  @author Florent Pruvost
#  @date 13-07-2012
#
###

cmake_minimum_required(VERSION 2.8)

# check if magma_dgetrf_incpiv_gpu is accessible in libmagma and activate it in chameleon
if ( CBLAS_FOUND AND LAPACKE_FOUND AND LAPACK_FOUND AND CUDA_FOUND AND CUDA_CUBLAS_LIBRARIES AND MAGMA_FOUND )
    set(CMAKE_REQUIRED_LIBRARIES "${CBLAS_LIBRARIES};${LAPACKE_LIBRARIES};${LAPACK_SEQ_LIBRARIES};${CUDA_LIBRARIES};${CUDA_CUBLAS_LIBRARIES};${MAGMA_LIBRARIES};${COREBLAS_LIBRARIES}")
    unset(MAGMA_DGETRF_INCPIV_GPU_FOUND CACHE)
    check_function_exists(magma_dgetrf_incpiv_gpu MAGMA_DGETRF_INCPIV_GPU_FOUND)
    if ( MAGMA_DGETRF_INCPIV_GPU_FOUND )
        message(STATUS "Set HAVE_MAGMA_GETRF_INCPIV_GPU")
        set(HAVE_MAGMA_GETRF_INCPIV_GPU 1)
    endif()
endif()

# Generate headers for all possible precisions
# --------------------------------------------
set(RUNTIME_HDRS_GENERATED "")
set(ZHDR
  include/runtime_codelet_z.h
)


# Define the list of headers
# --------------------------
set(RUNTIME_HDRS
    include/morse_starpu.h
    include/runtime_codelet_profile.h
    include/runtime_codelets.h
    include/runtime_profiling.h
    include/runtime_workspace.h
    ${RUNTIME_HDRS_GENERATED}
    )

# Force generation of headers
# ---------------------------
add_custom_target(runtime_starpu_include ALL SOURCES ${RUNTIME_HDRS})

# installation
# ------------
#install(FILES ${RUNTIME_HDRS}
#        DESTINATION include)

# Generate the morse common for all possible precisions
# -----------------------------------------------------
set(RUNTIME_COMMON_GENERATED "")
set(ZSRC
    control/runtime_zprofiling.c
    control/runtime_zlocality.c
    )

precisions_rules_py(RUNTIME_COMMON_GENERATED "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}"
                    TARGETDIR "control")

set(RUNTIME_COMMON
    control/runtime_async.c
    control/runtime_context.c
    control/runtime_control.c
    control/runtime_descriptor.c
    control/runtime_options.c
    control/runtime_profiling.c
    control/runtime_workspace.c
    codelets/codelet_dataflush.c
    ${RUNTIME_COMMON_GENERATED}
    )

set(flags_to_add "")
foreach(_prec ${CHAMELEON_PRECISION})
    set(flags_to_add "${flags_to_add} -DPRECISION_${_prec}")
endforeach()
set_source_files_properties(control/runtime_profiling.c PROPERTIES COMPILE_FLAGS "${flags_to_add}")

# Generate the morse sources for all possible precisions
# ------------------------------------------------------
set(RUNTIME_SRCS_GENERATED "")
set(ZSRC
    codelets/codelet_zcallback.c
    ${CODELETS_ZSRC}
    )

precisions_rules_py(RUNTIME_SRCS_GENERATED "${ZSRC}"
                    PRECISIONS "${CHAMELEON_PRECISION}"
                    TARGETDIR "codelets")

set(RUNTIME_SRCS
  ${RUNTIME_COMMON}
  ${RUNTIME_SRCS_GENERATED}
  )

# Add library
# -----------
add_library(chameleon_starpu ${RUNTIME_SRCS})
set_property(TARGET chameleon_starpu PROPERTY LINKER_LANGUAGE Fortran)
set_property(TARGET chameleon_starpu PROPERTY INSTALL_NAME_DIR "${CMAKE_INSTALL_PREFIX}/lib")

target_link_libraries(chameleon_starpu ${STARPU_LIBRARIES_DEP})
if (NOT CHAMELEON_SIMULATION)
    target_link_libraries(chameleon_starpu coreblas)
endif()
if(CHAMELEON_USE_CUDA AND NOT CHAMELEON_SIMULATION)
    target_link_libraries(chameleon_starpu cudablas)
endif(CHAMELEON_USE_CUDA AND NOT CHAMELEON_SIMULATION)

add_dependencies(chameleon_starpu
  chameleon_include
  coreblas_include
  control_include
  runtime_starpu_include
)
if (NOT CHAMELEON_SIMULATION)
    add_dependencies(chameleon_starpu coreblas_include)
endif()

if (CHAMELEON_USE_CUDA AND NOT CHAMELEON_SIMULATION)
    add_dependencies(chameleon_starpu cudablas_include)
endif()

# installation
# ------------
install(TARGETS chameleon_starpu
        DESTINATION lib)

###
### END CMakeLists.txt
###
