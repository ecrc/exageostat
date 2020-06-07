cmake_policy(SET CMP0057 NEW)

set( tests_files synthetic_dmle_test.c   synthetic_smle_test.c synthetic_sdmle_test.c
	 synthetic_dmle_locs_rwrapper_test.c synthetic_dmle_rwrapper_test.c real_csv_dmle_test.c
	 real_csv_smle_test.c real_csv_sdmle_test.c)

if(EXAGEOSTAT_USE_NETCDF)
	list(APPEND tests_files real_nc_dmle_1dlocs_test.c real_nc_dmle_2dlocs_test.c real_nc_sdmle_2dlocs_test.c real_nc_sdmle_1dlocs_test.c)
endif()

foreach(test_src ${tests_files})
    get_filename_component(test_exe ${test_src} NAME_WE)


if(EXAGEOSTAT_USE_CUDA)
#message(WARNING "USING CUDA")
	cuda_add_executable(${test_exe} ${test_src} examples.c)
        set_target_properties(${test_exe} PROPERTIES LINKER_LANGUAGE C )
else()
#message(WARNING "NOT USING CUDA")
	add_executable(${test_exe} ${test_src} examples.c)
endif()

    target_link_libraries(${test_exe} exageostat)
    #set_target_properties(test_${test_exe} PROPERTIES OUTPUT_NAME ${test_exe})
    list(APPEND exe_names ${test_exe})
endforeach()


#foreach(test_src ${tests_files})
#    get_filename_component(test_exe ${test_src} NAME_WE)
#    add_executable(${test_exe} ${test_src} examples.c)
#    target_link_libraries(${test_exe} exageostat)
#    #set_target_properties(test_${test_exe} PROPERTIES OUTPUT_NAME ${test_exe})
#    list(APPEND exe_names ${test_exe})
#endforeach()

file(DOWNLOAD "https://docs.google.com/uc?export=download&id=0B6mGsgjYG40kOGdNc2cyZzBjelk" "${CMAKE_CURRENT_BINARY_DIR}/marc-sample.tar.gz" EXPECTED_HASH MD5=2978769f72250eeff401f725158587a3)
file(DOWNLOAD "https://docs.google.com/uc?export=download&id=1Ib7W9Y5z7nqL2t9C8bbXFxL3J0E9f34U" "${CMAKE_CURRENT_BINARY_DIR}/wind-sample.tar.gz" EXPECTED_HASH MD5=2103ba786a3b07987ef70b4a37a94e9b)

add_custom_target( extract_tests_data ALL
                    DEPENDS ${exe_names}
                    COMMAND cmake -E tar zxf marc-sample.tar.gz
                    COMMAND cmake -E tar zxf wind-sample.tar.gz
  )
set( test_env  "STARPU_SILENT=1" )
set( test_env2 "STARPU_LIMIT_CPU_MEM=100" )
set( test_env3 "STARPU_DISK_SWAP=\"/tmp\"" )
set( test_env4 "STARPU_DISK_SWAP_BACKEND=unistd" )
set( test_env5 "STARPU_DISK_SWAP_SIZE=1024" )


#Test(1)
add_test( NAME synthetic_dmle_test_matern_stationary_ooc
          COMMAND numactl --interleave=all
           ./synthetic_dmle_test
          --test
          --N=6400
          --dts=512
          --ncores=${NUMBER_OF_CPU}
          --computation=exact
          --kernel=?:?:?
          --ikernel=1:0.1:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --opt_iters=10
          --opt_tol=4
          --ooc
          --kernel_fun=univariate_matern_stationary
)
set_tests_properties( synthetic_dmle_test_matern_stationary_ooc
                      PROPERTIES ENVIRONMENT "${test_env};${test_env2};${test_env3};${test_env4};${test_env5}" )

#Test(2)
add_test( NAME synthetic_dmle_test_ooc
          COMMAND numactl --interleave=all
           ./synthetic_dmle_test
          --test
          --N=6400
          --dts=512
          --ncores=${NUMBER_OF_CPU}
          --computation=exact
          --kernel=?:?:?
          --ikernel=1:0.1:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --opt_iters=10
          --opt_tol=4
          --ooc
          --kernel_fun=univariate_matern_stationary

)
set_tests_properties( synthetic_dmle_test_ooc
                      PROPERTIES ENVIRONMENT "${test_env};${test_env2};${test_env3};${test_env4};${test_env5}" )

#Test(3)
add_test( NAME synthetic_sdmle_test_ooc
          COMMAND numactl --interleave=all
           ./synthetic_sdmle_test
          --test
          --N=6400
          --dts=512
          --ncores=${NUMBER_OF_CPU}
          --computation=exact
    	  --diag_thick=1
          --kernel=?:?:?
          --ikernel=1:0.1:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --opt_iters=10
          --opt_tol=4
          --ooc
          --kernel_fun=univariate_matern_stationary

)
set_tests_properties( synthetic_sdmle_test_ooc
                      PROPERTIES ENVIRONMENT "${test_env};${test_env2};${test_env3};${test_env4};${test_env5}" )

#Test(4)
add_test( NAME synthetic_smle_test_ooc
          COMMAND numactl --interleave=all
           ./synthetic_smle_test
          --test
          --N=6400
          --dts=512
          --ncores=${NUMBER_OF_CPU}
          --computation=exact
          --kernel=?:?:?
          --ikernel=1:0.1:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --opt_iters=10
          --opt_tol=4
          --ooc
          --kernel_fun=univariate_matern_stationary
)
set_tests_properties( synthetic_smle_test_ooc
                      PROPERTIES ENVIRONMENT "${test_env};${test_env2};${test_env3};${test_env4};${test_env5}" )

#Test(5)

#add_test( NAME real_csv_dmle_test_t1
#          COMMAND numactl --interleave=all
#              ./real_csv_dmle_test
#              --dts=560
#              --ncores=${NUMBER_OF_CPU}
#              --computation=exact
#              --kernel=?:?:?
#              --olb=0.01:0.01:0.1
#              --oub=5:5:5
#              --obs_dir=./marc-sample/ppt.complete.Y001
#              --locs_file=./marc-sample/METAinfo
#              --predict=10
#              --kernel_fun=univariate_matern_stationary
#              --mspe
#)
#set_tests_properties( real_csv_dmle_test_t1
#                      PROPERTIES ENVIRONMENT "${test_env}")

#Test(6)
add_test( NAME synthetic_dmle_test_diag_approx_t1
          COMMAND numactl --interleave=all
          ./synthetic_dmle_test
          --test
          --N=19881
          --dts=960
          --ncores=${NUMBER_OF_CPU}
          --computation=diag_approx
          --diag_thick=5
          --kernel=?:?:?
          --ikernel=1:0.1:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --predict=10
          --lts=1900
          --acc=5
          --maxrank=300
          --opt_iters=5
          --opt_tol=2
          --predict=100
          --kernel_fun=univariate_matern_stationary
  )
set_tests_properties( synthetic_dmle_test_diag_approx_t1
                      PROPERTIES ENVIRONMENT "${test_env}")
#Test(7)
add_test( NAME synthetic_dmle_test_lr_approx_t1
          COMMAND numactl --interleave=all
          ./synthetic_dmle_test
          --test
          --N=8100
          --dts=960
          --ncores=${NUMBER_OF_CPU}
          --computation=lr_approx
          --kernel=?:?:?
          --ikernel=1:0.1:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --predict=100
          --lts=1000
          --acc=7
          --maxrank=400
          --opt_iters=2
          --opt_tol=4
          --kernel_fun=univariate_matern_stationary
  )
set_tests_properties( synthetic_dmle_test_lr_approx_t1
                      PROPERTIES ENVIRONMENT "${test_env}")
#Test(8)
#add_test( NAME synthetic_dmle_rwrapper_test_exact_t1
#         COMMAND numactl --interleave=all
#         ./synthetic_dmle_rwrapper_test
#         --test
#         --N=8100
#         --dts=560
#         --ncores=${NUMBER_OF_CPU}
#         --computation=exact
#         --kernel=?:?:?
#         --ikernel=1:0.03:0.5
#         --olb=0.01:0.01:0.01
#         --oub=5:5:5
#         --zvecs=1
#         --predict=10
#         --opt_iters=20
#         --opt_tol=4
#         --kernel_fun=univariate_matern_stationary
# )
#set_tests_properties( synthetic_dmle_rwrapper_test_exact_t1
#                      PROPERTIES ENVIRONMENT "${test_env}")

#Test(9)
#add_test( NAME synthetic_dmle_locs_rwrapper_test_exact_t1
#         COMMAND numactl --interleave=all
#         ./synthetic_dmle_locs_rwrapper_test
#         --dts=560
#          --ncores=${NUMBER_OF_CPU}
#         --computation=exact
#         --kernel=?:?:?
#         --ikernel=1:0.03:0.5
#         --olb=0.01:0.01:0.01
#         --oub=5:5:5
#         --zvecs=1
#         --locs_file=./marc-sample/METAinfo
#         --predict=10
#         --opt_iters=20
#         --opt_tol=4
#         --kernel_fun=univariate_matern_stationary
# )
#set_tests_properties( synthetic_dmle_locs_rwrapper_test_exact_t1
#                     PROPERTIES ENVIRONMENT "${test_env}")



if(EXAGEOSTAT_USE_NETCDF)
#Test(10)
#    add_test( NAME real_nc_dmle_1dlocs_test_exact
#             COMMAND numactl --interleave=all
#             ./real_nc_dmle_1dlocs_test
#              --dts=560
#             --ncores=${NUMBER_OF_CPU}
#             --computation=exact
#             --kernel=?:?:?
#             --olb=1.6:9.6:1.8
#             --oub=1.600004:9.625260:1.840517
#             --obs_dir=./20170901_UV_sample.nc
#             --locs_file=./20170901_UV_sample.nc
#             --predict=10
#             --dm=gc
#             --kernel_fun=univariate_matern_stationary
#     )
#set_tests_properties( real_nc_dmle_1dlocs_test_exact
#                     PROPERTIES ENVIRONMENT "${test_env};${test_env2};${test_env3};${test_env4};${test_env5}" )
#Test(11)
#    add_test( NAME real_nc_dmle_1dlocs_test_diag_approx
#             COMMAND numactl  --interleave=all
#             ./real_nc_dmle_1dlocs_test
#              --dts=560
#              --ncores=${NUMBER_OF_CPU}
#              --computation=diag_approx
#              --diag_thick=5
#              --kernel=?:?:?
#              --olb=1.6:9.6:1.8
#              --oub=1.600004:9.625260:1.840517
#              --obs_dir=./20170901_UV_sample.nc
#              --locs_file=./20170901_UV_sample.nc
#              --predict=10
#              --dm=gc
#              --kernel_fun=univariate_matern_stationary
#      )
#    set_tests_properties( real_nc_dmle_1dlocs_test_diag_approx
#                          PROPERTIES ENVIRONMENT "${test_env}")
endif()
#Test(12)
add_test( NAME synthetic_dmle_test_mloe_mmom
          COMMAND numactl --interleave=all
           ./synthetic_dmle_test
          --test
          --N=6400
          --dts=512
          --ncores=${NUMBER_OF_CPU}
          --computation=exact
          --kernel=?:?:?
          --ikernel=1:0.1:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --opt_iters=10
          --opt_tol=4
          --predict=10
    	  --mloe_mmom
          --kernel_fun=univariate_matern_stationary
)

set_tests_properties( synthetic_dmle_test_mloe_mmom
                      PROPERTIES ENVIRONMENT "${test_env}")

#Test(13)
add_test( NAME synthetic_dmle_test_bivariate_parsimonious
          COMMAND numactl --interleave=all

           ./synthetic_dmle_test
          --test
          --N=6400
          --dts=512
          --ncores=${NUMBER_OF_CPU}
          --computation=exact
          --kernel=?:?:?:?:?:?
          --ikernel=1:1:0.1:0.5::0.5:0.1
          --olb=0.01:0.01:0.01:0.01:0.01:0.01
          --oub=5:5:5:5:5:5
          --zvecs=1
          --opt_iters=10
          --opt_tol=4
          --predict=10
          --kernel_fun=bivariate_matern_parsimonious
)
set_tests_properties( synthetic_dmle_test_bivariate_parsimonious
                      PROPERTIES ENVIRONMENT "${test_env}")

#Test(14)
add_test( NAME synthetic_dmle_test_bivariate_parsimonious_tlr
          COMMAND numactl --interleave=all
           ./synthetic_dmle_test
          --test
          --N=1600
          --dts=320
          --ncores=${NUMBER_OF_CPU}
          --computation=lr_approx
          --kernel=?:?:?:?:?:?
          --ikernel=1:1:0.1:0.5::0.5:0.1
          --olb=0.01:0.01:0.01:0.01:0.01:0.01
          --oub=5:5:5:5:5:5
          --zvecs=1
          --opt_iters=5
          --opt_tol=4
		  --acc=5
		  --maxrank=400
		  --lts=600	
          --kernel_fun=bivariate_matern_parsimonious
          )
set_tests_properties( synthetic_dmle_test_bivariate_parsimonious_tlr
                      PROPERTIES ENVIRONMENT "${test_env}")


#Test(15)
add_test( NAME synthetic_dmle_test_bivariate_parsimonious_dst
          COMMAND numactl --interleave=all
           ./synthetic_dmle_test
          --test
          --N=1600
          --dts=320
          --ncores=${NUMBER_OF_CPU}
          --computation=diag_approx
          --kernel=?:?:?:?:?:?
          --ikernel=1:1:0.1:0.5::0.5:0.1
          --olb=0.01:0.01:0.01:0.01:0.01:0.01
          --oub=5:5:5:5:5:5

          --zvecs=1
          --opt_iters=10
          --opt_tol=4
          --diag_thick=2
          --kernel_fun=bivariate_matern_parsimonious
)
set_tests_properties( synthetic_dmle_test_bivariate_parsimonious_dst
                          PROPERTIES ENVIRONMENT "${test_env}")
