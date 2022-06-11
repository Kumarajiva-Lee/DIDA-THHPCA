project(gmcore LANGUAGES Fortran)
option(V_POLE "Stagger V on the two Poles" OFF)
option(R4 "Use single precision" OFF)
option(R16 "Use quadruple precision" OFF)
option(Detail_Time "Output Detila_time time_flow" OFF)
option(Detail_Calc_Time "Output Detila_Calc_time" OFF)
option(Detail_io_Time "Output Detila_io_time" OFF)
option(ENSURE_ORDER "Ensure computation order in parallel" ON)
option(Only_Gmcore  "Only Gmcore Run" OFF)
option(member_num_1 "Use member 1" OFF)
option(member_num_2 "Use member 2" OFF)
option(member_num_4 "Use member 4" OFF)
option(member_num_8 "Use member 8" ON)
option(member_num_16 "Use member 16" OFF)
option(member_num_32 "Use member 32" OFF)
option(member_num_64 "Use member 64" OFF)
option(member_num_128 "Use member 128" OFF)

if (V_POLE)
  add_definitions(-DV_POLE)
	set(CPPFLAGS ${CPPFLAGS} -DV_POLE)
	message(STATUS "Place v at Poles.")
endif ()
if (R4)
  add_definitions(-DREAL_KIND=4)
	set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=4)
  message(STATUS "Use single precision.")
elseif (R16)
  add_definitions(-DREAL_KIND=16)
	set(CPPFLAGS ${CPPFLAGS} -DREAL_KIND=16)
  message(STATUS "Use quadruple precision.")
else ()
  add_definitions(-DREAL_KIND=8)
endif ()
if (member_num_1)
  add_definitions(-Dmember_num=1)
	set(CPPFLAGS ${CPPFLAGS} -Dmember_num=1)
  message(STATUS "Use member 1")
elseif (member_num_2)
  add_definitions(-Dmember_num=2)
  set(CPPFLAGS ${CPPFLAGS} -Dmember_num=2)
  message(STATUS "Use member 2")
elseif (member_num_4)
  add_definitions(-Dmember_num=4)
  set(CPPFLAGS ${CPPFLAGS} -Dmember_num=4)
  message(STATUS "Use member 4")
elseif (member_num_8)
  add_definitions(-Dmember_num=8)
  set(CPPFLAGS ${CPPFLAGS} -Dmember_num=8)
  message(STATUS "Use member 8")
elseif (member_num_16)
  add_definitions(-Dmember_num=16)
  set(CPPFLAGS ${CPPFLAGS} -Dmember_num=16)
  message(STATUS "Use member 16")
elseif (member_num_32)
  add_definitions(-Dmember_num=32)
  set(CPPFLAGS ${CPPFLAGS} -Dmember_num=32)
  message(STATUS "Use member 32")
elseif (member_num_64)
  add_definitions(-Dmember_num=64)
  set(CPPFLAGS ${CPPFLAGS} -Dmember_num=64)
  message(STATUS "Use member 64")
elseif (member_num_128)
  add_definitions(-Dmember_num=128)
  set(CPPFLAGS ${CPPFLAGS} -Dmember_num=128)
  message(STATUS "Use member 128")
endif()
if (Detail_Time)
  add_definitions(-DDetail_Time)
  set(CPPFLAGS ${CPPFLAGS} -DDetail_Time)
  message(STATUS "Output Detila_time time_flow.")
endif ()
if (Detail_Calc_Time)
  add_definitions(-DDetail_Calc_Time)
  set(CPPFLAGS ${CPPFLAGS} -DDetail_Calc_Time)
  message(STATUS "Output Detila_Calc_time.")
endif ()
if (Detail_io_Time)
  add_definitions(-DDetail_io_Time)
  set(CPPFLAGS ${CPPFLAGS} -DDetail_io_Time)
  message(STATUS "Output Detila_io_time.")
endif ()
if (ENSURE_ORDER)
  add_definitions(-DENSURE_ORDER)
  set(CPPFLAGS ${CPPFLAGS} -DENSURE_ORDER)
  message(STATUS "Ensure computation order in parallel.")
endif ()
if (Only_Gmcore)
  add_definitions(-DOnly_Gmcore)
	set(CPPFLAGS ${CPPFLAGS} -DOnly_Gmcore)
	message(STATUS "Only Gmcore Run.")
endif ()


if(_is_sw_platform)
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -g -O3 -ffree-line-length-none ")
else()
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -ffpe-trap=invalid,zero,overflow")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O0 -ftrapuv -fpe0 -fcheck=all -g -fbacktrace")
    else ()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Ofast")
    endif ()
    # FIXME: OpenMP causes sigmentation faults.
    # set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fopenmp")
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -no-wrap-margin -O0 -check bounds -g -traceback -mkl")
    else ()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -traceback -no-wrap-margin -xCORE-AVX512  -O3 -align array64byte")
    endif ()
  endif ()
endif ()

if(_is_sw_platform)
  find_package(MPI)
  include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})
else()
  find_package(MPI)
  include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})
endif()

if(_is_sw_platform)
  if (DEFINED ENV{NETCDF_ROOT})
    include_directories("$ENV{NETCDF_ROOT}/include")
    link_directories("$ENV{NETCDF_ROOT}/lib")
  else ()
    message(FATAL_ERROR "Pkg-config could not find sw netcdf library!")
    message(FATAL_ERROR "Please set NETCDF_ROOT path or NETCDF path!")
  endif()
else()
  if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
    set(ENV{NETCDF_ROOT} $ENV{NETCDF})
  endif ()
  if (DEFINED ENV{NETCDF_ROOT})
    include_directories("$ENV{NETCDF_ROOT}/include")
    link_directories("$ENV{NETCDF_ROOT}/lib")
  else ()
    find_package(PkgConfig)
    if (PKG_CONFIG_FOUND)
      pkg_search_module(NETCDF REQUIRED netcdf)
      if (NETCDF_FOUND)
        include_directories(${NETCDF_INCLUDEDIR})
        link_directories(${NETCDF_LIBDIR})
      else ()
        message(FATAL_ERROR "Pkg-config could not find netcdf library!")
      endif ()
    else ()
      message(FATAL_ERROR "Unable to find pkg-config library!")
    endif ()
  endif ()
endif()

#add back~~!!!
#if (DEFINED ENV{MKLROOT})
#  add_definitions(-DHAS_MKL)
  # include_directories("$ENV{MKLROOT}/include")
  # link_directories("$ENV{MKLROOT}/lib")
#elseif (DEFINED ENV{MKL_INC} AND DEFINED ENV{MKL_LIB})
#  add_definitions(-DHAS_MKL)
  # include_directories("$ENV{MKL_INC}")
  # link_directories("$ENV{MKL_LIB}")
#endif()
#add back~~!!!


#if (DEFINED ENV{MKLROOT})
#  include_directories("$ENV{MKLROOT}/include")
#  link_directories("$ENV{MKLROOT}/lib")
#elseif (DEFINED ENV{MKL_INC} AND DEFINED ENV{MKL_LIB})
#  include_directories("$ENV{MKL_INC}")
#  link_directories("$ENV{MKL_LIB}")
#else ()
#  message (FATAL_ERROR "Unable to find MKL library!")
#endif ()



if(_is_sw_platform)
else()
  if (HAS_MKL)
    if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
      set(EXTERNAL_LIBS ${EXTERNAL_LIBS} lapack64)
    endif ()
  endif ()
endif()

#add_subdirectory(lib/container)
#add_subdirectory(lib/datetime)
#add_subdirectory(lib/string)
#add_subdirectory(lib/fiona)
#add_subdirectory(lib/flogger)

add_custom_command(
  OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/reduce_mod.F90
  COMMAND cpp ${CPPFLAGS} -P -C ${PROJECT_SOURCE_DIR}/core/reduce_mod.F90 > ${CMAKE_CURRENT_BINARY_DIR}/reduce_mod.F90
  MAIN_DEPENDENCY ${PROJECT_SOURCE_DIR}/core/reduce_mod.F90
  VERBATIM
)

set(source_files
  utils/const_mod.F90
  utils/formula_mod.F90
  utils/namelist_mod.F90
  utils/member_mod.F90
  utils/math_mod.F90
  utils/diag_state_mod.F90
  utils/performance_analysis_mod.F90
  utils/history_mod.F90
  utils/initial_mod.F90
  utils/restart_mod.F90
  utils/time_mod.F90
  utils/sphere_geometry_mod.F90
  utils/block_mod.F90
  utils/process_mod.F90
  utils/halo_mod.F90
  utils/allocator_mod.F90
  utils/parallel_types_mod.F90
  utils/parallel_mod.F90
  utils/tridiag_mkl_mod.F90
  utils/tridiag_spk_mod.F90
  utils/tridiag_mod.F90
  utils/debug_mod.F90
  core/time_schemes_mod.F90
  core/pgf/pgf_mod.F90
  core/pgf/pgf_swm_mod.F90
  core/pgf/pgf_lin97_mod.F90
  core/mesh_mod.F90
  core/vert_coord/sigma_coord_mod.F90
  core/vert_coord/hybrid_coord_test_mod.F90
  core/vert_coord/hybrid_coord_ecmwf_mod.F90
  core/vert_coord/hybrid_coord_mod.F90
  core/vert_coord/vert_coord_mod.F90
  core/refer_state/refer_state_types_mod.F90
  core/refer_state/refer_state_mod.F90
  core/refer_state/refer_state_wrf_mod.F90
  core/state_mod.F90
  core/static_mod.F90
  core/tend_mod.F90
  core/nh_mod.F90
  core/interp_mod.F90
  core/upwind_mod.F90
  core/weno_mod.F90
  core/operators_mod.F90
  core/filter_types_mod.F90
  core/filter_mod.F90
  core/damp/div_damp_mod.F90
  core/damp/vor_damp_mod.F90
  core/damp/smag_damp_mod.F90
  core/damp/damp_mod.F90
  core/damp/laplace_damp_mod.F90
  prepare/prepare_mod.F90
  prepare/topo_mod.F90
  prepare/latlon_interp_mod.F90
  prepare/vert_interp_mod.F90
  prepare/era5_reader_mod.F90
  prepare/mpas_reader_mod.F90
  prepare/bkg_mod.F90
  prepare/waccm_reader_mod.F90
  prepare/openmars_reader_mod.F90
  core/gmcore_mod.F90
  physics/test_forcing_mod.F90
  tests/swm/rossby_haurwitz_wave_test_mod.F90
  tests/swm/mountain_zonal_flow_test_mod.F90
  tests/swm/jet_zonal_flow_test_mod.F90
  tests/swm/steady_geostrophic_flow_test_mod.F90
  tests/swm/cross_pole_flow_test_mod.F90
  tests/swm/shallow_water_waves_test_mod.F90
  tests/swm/vortex_erosion_test_mod.F90
  tests/steady_state_test_mod.F90
  tests/rossby_haurwitz_wave_3d_test_mod.F90
  tests/mountain_wave_test_mod.F90
  tests/baroclinic_wave_test_mod.F90
  tests/held_suarez_test_mod.F90
  tests/steady_state_pgf_test_mod.F90
  tests/ksp15_test_mod.F90
  tests/dcmip31_test_mod.F90
  tests/mars_cold_run_mod.F90
  ../LIB/quadpack.f90
  gmcore_mct/gmcore_coupler.F90
)

add_library(atm ${source_files})
#target_link_libraries(atm netcdff coupler fortran_container gptl fiona flogger fortran_datetime ${MPI_Fortran_LIBRARIES} ${MKL_LIBRARIES})
target_link_libraries(atm netcdff coupler fortran_container gptl fiona flogger fortran_datetime mkl_intel_lp64 mkl_sequential mkl_core)

#add_executable(gmcore_swm_driver.exe  drivers/gmcore_swm_driver.F90)
#target_link_libraries(gmcore_swm_driver.exe gmcore)

#add_executable(gmcore_driver.exe  drivers/gmcore_driver.F90)
#target_link_libraries(gmcore_driver.exe gmcore)
