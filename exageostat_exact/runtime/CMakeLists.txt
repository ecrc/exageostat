# List of codelets required by all runtimes
# -----------------------------------------
set(CODELETS_DSRC
	codelets/codelet_dcmg.c
	codelets/codelet_ddotp.c
	codelets/codelet_dgemv.c
	codelets/codelet_dmdet.c
	codelets/codelet_dmse.c
	codelets/codelet_dzcpy.c
    )

# Check for the subdirectories
# ----------------------------
if( CHAMELEON_SCHED_QUARK )
    add_subdirectory(quark)
elseif( CHAMELEON_SCHED_PARSEC )
    add_subdirectory(parsec)
elseif( CHAMELEON_SCHED_STARPU )
    add_subdirectory(starpu)
endif()

###
### END CMakeLists.txt
###
