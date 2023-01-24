cmake_policy(SET CMP0057 NEW)

set(tests_files synthetic_dmle_test.c synthetic_smle_test.c synthetic_sdmle_test.c
	synthetic_dmle_locs_rwrapper_test.c synthetic_dmle_rwrapper_test.c real_csv_dmle_test.c
	synthetic_sdmle_rwrapper_test.c real_csv_smle_test.c real_csv_sdmle_test.c)

if (EXAGEOSTAT_USE_NETCDF)
	list(APPEND tests_files real_nc_dmle_1dlocs_test.c real_nc_dmle_2dlocs_test.c real_nc_sdmle_2dlocs_test.c real_nc_sdmle_1dlocs_test.c)
endif ()

foreach (test_src ${tests_files})
	get_filename_component(test_exe ${test_src} NAME_WE)
	if (EXAGEOSTAT_USE_CUDA)
		cuda_add_executable(${test_exe} ${test_src} examples.c)
		set_target_properties(${test_exe} PROPERTIES LINKER_LANGUAGE C)
	else ()
		add_executable(${test_exe} ${test_src} examples.c)
	endif ()

	target_link_libraries(${test_exe} exageostat)
	list(APPEND exe_names ${test_exe})
endforeach ()

file(DOWNLOAD "https://docs.google.com/uc?export=download&id=13B0kNFGdTNV3LQOsutHCxtQTSoxTiz1S" "${CMAKE_CURRENT_BINARY_DIR}/marc-sample.tar.gz")
file(DOWNLOAD "https://docs.google.com/uc?export=download&id=1Ib7W9Y5z7nqL2t9C8bbXFxL3J0E9f34U" "${CMAKE_CURRENT_BINARY_DIR}/wind-sample.tar.gz")

add_custom_target(extract_tests_data ALL
	DEPENDS ${exe_names}
	COMMAND cmake -E tar zxf marc-sample.tar.gz
	COMMAND cmake -E tar zxf wind-sample.tar.gz
	)
set(test_env "STARPU_SILENT=1")
set(test_env2 "STARPU_LIMIT_CPU_MEM=100")
set(test_env3 "STARPU_DISK_SWAP=\"/tmp\"")
set(test_env4 "STARPU_DISK_SWAP_BACKEND=unistd")
set(test_env5 "STARPU_DISK_SWAP_SIZE=1024")

#Test(1): Univariate Matern Gaussian stationary kernel in dense -- double precision
add_test(NAME synthetic_dmle_test_matern_stationary
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=6400
	--dts=512
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?
	--ikernel=1:0.1:0.5
	--olb=0.1:0.1:0.1
	--oub=5:5:5
	--zvecs=1
	--opt_iters=10
	--opt_tol=4
	--kernel_fun=univariate_matern_stationary
	)
set_tests_properties(synthetic_dmle_test_matern_stationary
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(2): Univariate Matern Gaussian stationary kernel in dense -- single precision
add_test(NAME synthetic_smle_test_matern_stationary
	COMMAND numactl --interleave=all
	./synthetic_smle_test
	--test
	--N=6400
	--dts=512
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?
	--ikernel=1:0.1:0.5
	--olb=0.1:0.1:0.1
	--oub=5:5:5
	--zvecs=1
	--opt_iters=5
	--opt_tol=4
	--kernel_fun=univariate_matern_stationary
	)
set_tests_properties(synthetic_smle_test_matern_stationary
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(3): Univariate Matern Gaussian stationary kernel with Nugget in dense -- double precision
add_test(NAME synthetic_dmle_test_matern_stationary_nugget
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=6400
	--dts=512
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?:?
	--ikernel=1:0.1:0.5:0.1
	--olb=0.1:0.1:0.1:0.01
	--oub=5:5:5:5
	--zvecs=1
	--opt_iters=10
	--opt_tol=4
	--kernel_fun=univariate_matern_nuggets_stationary

	)
set_tests_properties(synthetic_dmle_test_matern_stationary_nugget
	PROPERTIES ENVIRONMENT "${test_env}")


##Test(4): Univariate Matern Gaussian stationary kernel with Nugget in dense -- single precision
#add_test(NAME synthetic_smle_test_matern_stationary_nugget
#	COMMAND numactl --interleave=all
#	./synthetic_smle_test
#	--test
#	--N=6400
#	--dts=512
#	--ncores=${NUMBER_OF_CPU}
#	--computation=exact
#	--kernel=?:?:?:?
#	--ikernel=1:0.1:0.5:0.1
#	--olb=0.01:0.01:0.01:0.01
#	--oub=5:5:5:5
#	--zvecs=1
#	--opt_iters=10
#	--opt_tol=4
#	--kernel_fun=univariate_matern_nuggets_stationary


#	)
#set_tests_properties(synthetic_smle_test_matern_stationary_nugget
#	PROPERTIES ENVIRONMENT "${test_env}")

#Test(5): Univariate Matern Gaussian stationary kernel with Nugget in diag_approx -- double precision
add_test(NAME synthetic_dmle_test_matern_stationary_nugget_diag_approx
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=6400
	--dts=512
	--ncores=${NUMBER_OF_CPU}
	--computation=diag_approx
	--kernel=?:?:?:?
	--ikernel=1:0.1:0.5:0.1
	--olb=0.01:0.01:0.01:0.01
	--oub=5:5:5:5
	--zvecs=1
	--opt_iters=10
	--opt_tol=4
	--diag_thick=1
	--kernel_fun=univariate_matern_nuggets_stationary

	)
set_tests_properties(synthetic_dmle_test_matern_stationary_nugget_diag_approx
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(6): Matern Gaussian stationary kernel in dense -- mixed-precision
#add_test(NAME synthetic_sdmle_test
#        COMMAND numactl --interleave=all
#        ./synthetic_sdmle_test
#        --test
#        --N=6400
#        --dts=512
#        --ncores=${NUMBER_OF_CPU}
#        --computation=exact
#        --diag_thick=1
#        --kernel=?:?:?
#        --ikernel=1:0.1:0.5
#        --olb=0.1:0.1:0.1
#        --oub=5:5:5
#        --zvecs=1
#        --opt_iters=10
#        --opt_tol=4
#        --kernel_fun=univariate_matern_stationary
#        )
#
#set_tests_properties(synthetic_sdmle_test
#        PROPERTIES ENVIRONMENT "${test_env}")#;${test_env2};${test_env3};${test_env4};${test_env5}" )

#Test(7): Univariate Matern Gaussian stationary kernel in dense -- double precision
add_test(NAME synthetic_dmle_test_matern_stationary_ooc
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=6400
	--dts=512
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?
	--ikernel=1:0.1:0.5
	--olb=0.1:0.1:0.1
	--oub=5:5:5
	--zvecs=1
	--opt_iters=10
	--opt_tol=4
	--kernel_fun=univariate_matern_stationary
	--ooc
	)
set_tests_properties(synthetic_dmle_test_matern_stationary_ooc
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(8): Univariate Matern Gaussian stationary kernel in diag_approx -- double precision
add_test(NAME synthetic_dmle_test_diag_approx_t1
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=8100
	--dts=512
	--ncores=${NUMBER_OF_CPU}
	--computation=diag_approx
	--diag_thick=2
	--kernel=?:?:?
	--ikernel=1:0.1:0.5
	--olb=0.01:0.01:0.01
	--oub=3:1:3
	--zvecs=1
	--predict=10
	--opt_iters=5
	--opt_tol=2
	--kernel_fun=univariate_matern_stationary

	)
set_tests_properties(synthetic_dmle_test_diag_approx_t1
	PROPERTIES ENVIRONMENT "${test_env}")

if (EXAGEOSTAT_USE_HICMA)
	#Test(9): Univariate Matern Gaussian stationary kernel in lr_approx -- double precision
	add_test(NAME synthetic_dmle_test_lr_approx_t1
		COMMAND numactl --interleave=all
		./synthetic_dmle_test
		--test
		--N=1600
		--dts=320
		--ncores=${NUMBER_OF_CPU}
		--computation=lr_approx
		--kernel=1:?:0.5
		--ikernel=1:0.1:0.5
		--olb=0.1:0.1:0.1
		--oub=5:5:5
		--zvecs=1
		--predict=100
		--opt_tol=4
		--mspe
		--lts=1000
		--acc=7
		--maxrank=500
		--opt_iters=2
		--opt_tol=4
		--kernel_fun=univariate_matern_stationary
		)
	set_tests_properties(synthetic_dmle_test_lr_approx_t1
		PROPERTIES ENVIRONMENT "${test_env}")
endif ()

#Test(10): Univarate Matern Gaussian stationary kernel in dense + MLOE and MMOM -- double precision
add_test(NAME synthetic_dmle_test_mloe_mmom
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=6400
	--dts=512
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?
	--ikernel=1:0.1:0.5
	--olb=0.1:0.1:0.1
	--oub=5:5:5
	--zvecs=1
	--opt_iters=10
	--opt_tol=4
	--predict=10
	--mloe_mmom
	--kernel_fun=univariate_matern_stationary
	)
set_tests_properties(synthetic_dmle_test_mloe_mmom
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(11): Bivariate Parsimonious Matern Gaussian stationary kernel in dense -- double precision
add_test(NAME synthetic_dmle_test_bivariate_parsimonious
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=6400
	--dts=512
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?:?:?:?
	--ikernel=1:1:0.1:0.5:0.5:0.1
	--olb=0.1:0.1:0.1:0.1:0.1:0.1
	--oub=5:5:5:5:5:5
	--zvecs=1
	--opt_iters=5
	--opt_tol=4
	--predict=10
	--mspe
	--kernel_fun=bivariate_matern_parsimonious
	)
set_tests_properties(synthetic_dmle_test_bivariate_parsimonious
	PROPERTIES ENVIRONMENT "${test_env}")

if (EXAGEOSTAT_USE_HICMA)
	#Test(12): Bivariate Parsimonious Matern Gaussian stationary kernel in lr_approx -- double precision
	#    add_test(NAME synthetic_dmle_test_bivariate_parsimonious_tlr
	#            COMMAND numactl --interleave=all
	#            ./synthetic_dmle_test
	#            --test
	#            --N=1600
	#            --dts=320
	#            --ncores=${NUMBER_OF_CPU}
	#            --computation=lr_approx
	#            --kernel=?:?:?:?:?:?
	#            --ikernel=1:1:0.1:0.5:0.5:0.1
	#            --olb=0.1:0.1:0.1:0.1:0.1:0.1
	#            --oub=5:5:5:5:5:5
	#            --zvecs=1
	#            --opt_iters=5
	#            --opt_tol=4
	#            --acc=5
	#            --maxrank=400
	#            --lts=600
	#            --kernel_fun=bivariate_matern_parsimonious
	#            )
	#    set_tests_properties(synthetic_dmle_test_bivariate_parsimonious_tlr
	#            PROPERTIES ENVIRONMENT "${test_env}")
endif ()
#Test(13): Bivariate Parsimonious Matern Gaussian stationary kernel in diag_approx -- double precision
add_test(NAME synthetic_dmle_test_bivariate_parsimonious_dst
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=1600
	--dts=320
	--ncores=${NUMBER_OF_CPU}
	--computation=diag_approx
	--kernel=?:?:?:?:?:?
	--ikernel=1:1:0.1:0.5:0.5:0.1
	--olb=0.01:0.01:0.01:0.01:0.01:0.01
	--oub=2:2:2:2:2:2
	--zvecs=1
	--opt_iters=7
	--opt_tol=4
	--diag_thick=2
	--kernel_fun=bivariate_matern_parsimonious
	)
set_tests_properties(synthetic_dmle_test_bivariate_parsimonious_dst
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(14): Bivariate Parsimonious Matern Gaussian stationary kernel in diag_approx + MLOE and MMOM -- double precision
add_test(NAME synthetic_dmle_test_bivariate_parsimonious_mloe_mmom
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=1600
	--dts=320
	--ncores=${NUMBER_OF_CPU}
	--computation=diag_approx
	--kernel=?:?:?:?:?:?
	--ikernel=1:1:0.1:0.5:0.5:0.1
	--olb=0.1:0.1:0.1:0.1:0.1:0.1
	--oub=5:5:5:5:5:5
	--zvecs=1
	--opt_iters=4
	--opt_tol=4
	--kernel_fun=bivariate_matern_parsimonious
	--predict=3
	--diag_thick=1
	--mloe_mmom
	)
set_tests_properties(synthetic_dmle_test_bivariate_parsimonious_mloe_mmom
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(15): Bivariate Flexible Matern Gaussian stationary kernel in dense -- double precision
add_test(NAME synthetic_dmle_test_bivariate_flexible
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=8100
	--dts=320
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?:?:?:?:?:?:?:?:?
	--ikernel=0.3:0.6:0.01:0.3:0.9:0.9:0.05:0.3:1.5:0.9:0.99
	--olb=0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01
	--oub=50:50:50:50:50:50:50:50:50:50:50
	--zvecs=1
	--opt_iters=4
	--opt_tol=4
	--kernel_fun=bivariate_matern_flexible
	)
set_tests_properties(synthetic_dmle_test_bivariate_flexible
	PROPERTIES ENVIRONMENT "${test_env}")


#Test(16): Bivariate Flexible Matern Gaussian stationary kernel in diag_approx -- double precision
add_test(NAME synthetic_dmle_test_bivariate_flexible_diag_approx
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=8100
	--dts=320
	--ncores=${NUMBER_OF_CPU}
	--computation=diag_approx
	--kernel=?:?:?:?:?:?:?:?:?:?:?
	--ikernel=0.3:0.6:0.01:0.3:0.9:0.9:0.05:0.3:1.5:0.9:0.99
	--olb=0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01
	--oub=50:50:50:50:50:50:50:50:50:50:50
	--zvecs=1
	--opt_iters=4
	--opt_tol=4
	--diag_thick=2
	--kernel_fun=bivariate_matern_flexible
	)
set_tests_properties(synthetic_dmle_test_bivariate_flexible_diag_approx
	PROPERTIES ENVIRONMENT "${test_env}")


#Test(17): Univariate space-time Matern Gaussian stationary kernel in exact -- double precision
add_test(NAME synthetic_dmle_test_space_time
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=1600
	--dts=320
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?:?:?:?:0
	--ikernel=1:1:0.1:0.5:0.5:0.1:0
	--olb=0.1:0.1:0.1:0.1:0.1:0.1:0
	--oub=5:5:5:5:5:5:0
	--zvecs=1
	--opt_iters=4
	--opt_tol=20
	--kernel_fun=univariate_spacetime_matern_stationary
	--dim=st
	--time_slots=5
	)
set_tests_properties(synthetic_dmle_test_space_time
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(18): Trivariate Parsimonious Matern Gaussian stationary kernel in exact -- double precision
add_test(NAME synthetic_dmle_test_trivariate_matern_parsimonious
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=1600
	--dts=300
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?:?:?:?:?:?:?:?
	--ikernel=1:1:1:0.1:0.5:1:1.5:0.1:0.1:0
	--olb=0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1
	--oub=5:5:5:5:1:1:5:5:5:5
	--zvecs=1
	--opt_iters=8
	--opt_tol=5
	--diag_thick=2
	--kernel_fun=trivariate_matern_parsimonious
	)
set_tests_properties(synthetic_dmle_test_bivariate_parsimonious
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(19): Trivariate Parsimonious Matern Gaussian stationary kernel in diag_approx -- double precision
add_test(NAME synthetic_dmle_test_trivariate_matern_parsimonious_diag_approx
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=1600
	--dts=300
	--ncores=${NUMBER_OF_CPU}
	--computation=diag_approx
	--kernel=?:?:?:?:?:?:?:?:?:?
	--ikernel=1:1:1:0.1:0.5:1:1.5:0.1:0.1:0
	--olb=0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1
	--oub=5:5:5:5:1:1:5:5:5:5
	--zvecs=1
	--opt_iters=4
	--opt_tol=5
	--diag_thick=2
	--kernel_fun=trivariate_matern_parsimonious
	)
set_tests_properties(synthetic_dmle_test_trivariate_matern_parsimonious_diag_approx
	PROPERTIES ENVIRONMENT "${test_env}")


#Test(20): Trivariate Parsimonious Matern Gaussian stationary kernel in dense (optimization with profiling) -- double precision
add_test(NAME synthetic_dmle_test_trivariate_matern_parsimonious_profile
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=8100
	--dts=300
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?:?:?:?:?:?:?:?
	--ikernel=1:1:1:0.1:0.5:1:1.5:0.1:0.1:0
	--olb=0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1:0.1
	--oub=5:5:5:5:1:1:5:5:5:5
	--zvecs=1
	--opt_iters=5
	--opt_tol=5
	--kernel_fun=trivariate_matern_parsimonious_profile
	)
set_tests_properties(synthetic_dmle_test_trivariate_matern_parsimonious_profile
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(21): Univariate Gaussian non-stationary kernel in dense -- double precision
add_test(NAME synthetic_dmle_test_univariate_matern_non_stat
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=6400
	--dts=360
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?:?:?:?:?:?:?
	--ikernel=0.04:1.57:0.33:-1:0.8:0.1:-0.5:0.5
	--olb=0.01:0.01:0.01:0.01:0.01:0.01:0.01:0.01
	--oub=5:5:5:5:5:5:5:5
	--zvecs=1
	--opt_iters=1
	--opt_tol=4
	--kernel_fun=univariate_matern_non_stat
	)
set_tests_properties(synthetic_dmle_test_univariate_matern_non_stat
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(22): Univariate non-Gaussian stationary kernel in dense -- double precision
add_test(NAME synthetic_dmle_test_univariate_matern_non_gaussian
	COMMAND numactl --interleave=all
	./synthetic_dmle_test
	--test
	--N=6400
	--dts=320
	--ncores=${NUMBER_OF_CPU}
	--computation=exact
	--kernel=?:?:?:?:?:?
	--ikernel=7.0711:1:0:2:0:0
	--olb=0.01:0.01:-5:0.1:-2:0
	--oub=15:5:5:5:2:2
	--zvecs=1
	--opt_iters=1
	--opt_tol=5
	--kernel_fun=univariate_matern_non_gaussian
	)
set_tests_properties(synthetic_dmle_test_univariate_matern_non_gaussian
	PROPERTIES ENVIRONMENT "${test_env}")

#Test(23): Univariate non-Gaussian stationary kernel in diag_approx -- double precision
#add_test(NAME synthetic_dmle_test_univariate_matern_non_gaussian_diag_approx
#        COMMAND numactl --interleave=all
#        ./synthetic_dmle_test
#        --test
#        --N=6400
#        --dts=320
#        --ncores=${NUMBER_OF_CPU}
#        --computation=diag_approx
#        --kernel=?:?:?:?:?:?
#        --ikernel=7.0711:1:0:2:0:0
#        --olb=0.01:0.01:-5:0.1:-2:0
#        --oub=15:5:5:5:2:2
#        --zvecs=1
#        --opt_iters=1
#       --opt_tol=5
#        --diag_thick=2
#        --kernel_fun=univariate_matern_non_gaussian
#        )
#set_tests_properties(synthetic_dmle_test_univariate_matern_non_gaussian_diag_approx
#        PROPERTIES ENVIRONMENT "${test_env}")

if (EXAGEOSTAT_USE_HICMA)
	#Test(24): Univariate non-Gaussian stationary kernel in lr_approx -- double precision
	add_test(NAME synthetic_dmle_test_univariate_matern_non_gaussian_lr_approx
		COMMAND numactl --interleave=all
		./synthetic_dmle_test
		--test
		--N=6400
		--dts=320
		--ncores=1
		--computation=lr_approx
		--kernel=7.0711:1:0:2:0:0
		--ikernel=7.0711:1:0:2:0:0
		--olb=0.01:0.01:-5:0.1:-2:0
		--oub=15:5:5:5:2:2
		--zvecs=1
		--opt_iters=1
		--opt_tol=5
		--acc=7
		--lts=512
		--maxrank=200
		--kernel_fun=univariate_matern_non_gaussian
		)
	set_tests_properties(synthetic_dmle_test_univariate_matern_non_gaussian_lr_approx
		PROPERTIES ENVIRONMENT "${test_env}")
endif ()
