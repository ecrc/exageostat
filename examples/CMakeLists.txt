cmake_policy(SET CMP0057 NEW)

set( tests_files zgen_mle_test.c zgen_mle_testr.c mle_flatfile_real.c zgen_mle_givenlocs_testr.c)

if(EXAGEOSTAT_USE_NETCDF)
	list(APPEND tests_files mle_nc_real_1dlocs.c mle_nc_real_2dlocs.c)
endif()

foreach(test_src ${tests_files})
    get_filename_component(test_exe ${test_src} NAME_WE)
    add_executable(${test_exe} ${test_src})
    target_link_libraries(${test_exe} exageostat)
    #set_target_properties(test_${test_exe} PROPERTIES OUTPUT_NAME ${test_exe})
    list(APPEND exe_names ${test_exe})
endforeach()

file(DOWNLOAD "https://docs.google.com/uc?export=download&id=0B6mGsgjYG40kOGdNc2cyZzBjelk" "${CMAKE_CURRENT_BINARY_DIR}/marc-sample.tar.gz" EXPECTED_HASH MD5=2978769f72250eeff401f725158587a3)
file(DOWNLOAD "https://docs.google.com/uc?export=download&id=1Ib7W9Y5z7nqL2t9C8bbXFxL3J0E9f34U" "${CMAKE_CURRENT_BINARY_DIR}/wind-sample.tar.gz" EXPECTED_HASH MD5=2103ba786a3b07987ef70b4a37a94e9b)

add_custom_target( extract_tests_data ALL
                    DEPENDS ${exe_names}
                    COMMAND cmake -E tar zxf marc-sample.tar.gz
                    COMMAND cmake -E tar zxf wind-sample.tar.gz
  )
set( test_env "STARPU_SILENT=1" )
set( test_env2 "STARPU_LIMIT_CPU_MEM=100" )
set( test_env3 "STARPU_DISK_SWAP=\"/tmp\"" )
set( test_env4 "STARPU_DISK_SWAP_BACKEND=unistd" )
set( test_env5  "STARPU_DISK_SWAP_SIZE=1024" )

add_test( NAME zgen_mle_test_ooc
          COMMAND numactl --interleave=all
           ./zgen_mle_test
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
)
set_tests_properties( zgen_mle_test_ooc
                      PROPERTIES ENVIRONMENT "${test_env};${test_env2};${test_env3};${test_env4};${test_env5}" )


add_test( NAME mle_flatfile_real_t1
          COMMAND numactl --interleave=all
              ./mle_flatfile_real
              --dts=560
              --ncores=${NUMBER_OF_CPU}
              --computation=exact
              --kernel=?:?:?
              --olb=0.01:0.01:0.1
              --oub=5:5:5
              --obs_dir=./marc-sample/ppt.complete.Y001
              --locs_file=./marc-sample/METAinfo
              --predict=10
)
set_tests_properties( mle_flatfile_real_t1
                      PROPERTIES ENVIRONMENT "${test_env}")

add_test( NAME zgen_mle_test_diag_approx_t1
          COMMAND numactl --interleave=all
          ./zgen_mle_test
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
          --opt_iters=10
          --opt_tol=2
          --predict=100
  )
set_tests_properties( zgen_mle_test_diag_approx_t1
                      PROPERTIES ENVIRONMENT "${test_env}")

add_test( NAME zgen_mle_test_lr_approx_t1
          COMMAND numactl --interleave=all
          ./zgen_mle_test
          --test
          --N=19881
          --dts=960
          --ncores=${NUMBER_OF_CPU}
          --computation=lr_approx
          --kernel=?:?:?
          --ikernel=1:0.1:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --predict=10
          --lts=1900
          --acc=7
          --maxrank=400
          --opt_iters=10
          --opt_tol=4
  )
set_tests_properties( zgen_mle_test_lr_approx_t1
                      PROPERTIES ENVIRONMENT "${test_env}")

add_test( NAME zgen_mle_testr_exact_t1
          COMMAND numactl --interleave=all
          ./zgen_mle_testr
          --test
          --N=8100
          --dts=560
          --ncores=${NUMBER_OF_CPU}
          --computation=exact
          --kernel=?:?:?
          --ikernel=1:0.03:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --predict=10
          --opt_iters=20
          --opt_tol=4
  )
set_tests_properties( zgen_mle_testr_exact_t1
                      PROPERTIES ENVIRONMENT "${test_env}")


add_test( NAME zgen_mle_givenlocs_testr_exact_t1
          COMMAND numactl --interleave=all
          ./zgen_mle_givenlocs_testr
          --dts=560
          --ncores=${NUMBER_OF_CPU}
          --computation=exact
          --kernel=?:?:?
          --ikernel=1:0.03:0.5
          --olb=0.01:0.01:0.01
          --oub=5:5:5
          --zvecs=1
          --locs_file=./marc-sample/METAinfo
          --predict=10
          --opt_iters=20
          --opt_tol=4
  )
set_tests_properties( zgen_mle_givenlocs_testr_exact_t1
                      PROPERTIES ENVIRONMENT "${test_env}")


if(EXAGEOSTAT_USE_NETCDF)
    add_test( NAME mle_nc_real_1dlocs_exact
              COMMAND numactl --interleave=all
              ./mle_nc_real_1dlocs
              --dts=560
              --ncores=${NUMBER_OF_CPU}
              --computation=exact
              --kernel=?:?:?
              --olb=1.6:9.6:1.8
              --oub=1.600004:9.625260:1.840517
              --obs_dir=./20170901_UV_sample.nc
              --locs_file=./20170901_UV_sample.nc
              --predict=10
              --dm=gc
      )
    add_test( NAME mle_nc_real_1dlocs_diag_approx
              COMMAND numactl  --interleave=all
              ./mle_nc_real_1dlocs
              --dts=560
              --ncores=${NUMBER_OF_CPU}
              --computation=diag_approx
              --diag_thick=5
              --kernel=?:?:?
              --olb=1.6:9.6:1.8
              --oub=1.600004:9.625260:1.840517
              --obs_dir=./20170901_UV_sample.nc
              --locs_file=./20170901_UV_sample.nc
              --predict=10
              --dm=gc
      )
    set( test_env "STARPU_SILENT=1" )
    set_tests_properties( mle_nc_real_1dlocs_exact
                          PROPERTIES ENVIRONMENT "${test_env}")
    set_tests_properties( mle_nc_real_1dlocs_diag_approx
                          PROPERTIES ENVIRONMENT "${test_env}")
endif()
