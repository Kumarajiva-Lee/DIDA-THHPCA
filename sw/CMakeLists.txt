cmake_minimum_required(VERSION 3.0)

project (dida.exe C Fortran)

MESSAGE("Current source dir is" ${STACCATO_SOURCE_DIR})
#------------------------------------------------------------------------------------#
# 1. Options
# 2. SET
# 3. Compiler
# 4. MPI
# 5. NETCDF
# 6. MKL
# 7. include
# 8. Subdirectories
# 9. Target_link
#10. UNITTEST
#------------------------------------------------------------------------------------#

#####################################################################################
#  1. Options
#####################################################################################
set(CMAKE_VERBOSE_MAKEFILE ON)

option(SW "application run on sw machine" ON)
option(MASTER "xmath, libm run on sw master" OFF)
option(AUTOGPTL "use gptl auto instrumentation analysis" OFF)

option(LOG_MINFO "atm unittest on" ON)
option(RSL_LOG "each rank outputs to its own file" OFF)

option(ATM_COUPLING "atm coupling on" ON)
option(DA_COUPLING "atm coupling on" ON)

option(ATM_TEST "atm unittest on" ON)

option(DAOUT_NC "da output nctcdf" OFF)

option(LON_SLICE "redis use slice lon" ON)
option(REDIS_PIPE_W "redis use slice lon" OFF)
option(REDIS_PIPE_R "redis use slice lon" OFF)

option(fiona_TEST "fiona unittest" OFF)

#####################################################################################
#  2. SET
#####################################################################################
if(SW)
  set(_is_sw_platform 1)
else()
  set(_is_sw_platform 0)
endif()

if(MASTER)
  set(_is_sw_master 1)
else()
  set(_is_sw_master 0)
endif()

if(LOG_MINFO)
  add_definitions(-DLOG_MINFO=1)
endif()

if(RSL_LOG)
  add_definitions(-DRSL_LOG=1)
endif()

if(ATM_COUPLING)
  add_definitions(-DATM_COUPLING=1)
endif()

if(DA_COUPLING)
  add_definitions(-DDA_COUPLING=1)
endif()

if(DAOUT_NC)
  add_definitions(-DDAOUT_NC=1)
endif()

if(LON_SLICE)
  add_definitions(-DLON_SLICE=1)
endif()

if(REDIS_PIPE_W)
    add_definitions(-DREDIS_PIPE_W=1)
endif()

if(REDIS_PIPE_R)
    add_definitions(-DREDIS_PIPE_R=1)
endif()

####################################################################################
# 3. Compiler
####################################################################################
if(_is_sw_platform)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -g -O3 -mieee")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -std=c++11 -msimd")
else()
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3 ")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -std=c++11 ")
endif()

if(_is_sw_platform)
  #SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -finit-local-zero -mftz -g -O2 -ffree-line-length-none -fbacktrace -fcheck=all -mieee  -ffpe-trap=invalid,zero,overflow")
  SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -finit-local-zero -mftz -g -O3 -ffree-line-length-none -fcheck=no-array-temps -mieee ")
  SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /home/export/online1/mdt00/shisuan/swthhpca/liyiyuan/dida-sw/ATM/swdsl/slave/prepare_slave.o  /home/export/online1/mdt00/shisuan/swthhpca/liyiyuan/dida-sw/ATM/swdsl/host/c_host.o -lm_slave")
  #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L/usr/sw/lib -lmpiP -lbfd -liberty -lz")

  if(AUTOGPTL)
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -rdynamic -finstrument-functions -openmp")
    #SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -L../LIB/gptl/lib -lgptl -lunwind -lunwind-ptrace  -lunwind-setjmp -lunwind-generic -Wl,--eh-frame-hdr ")
    SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--eh-frame-hdr ")
  endif()
else()
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    if (CMAKE_BUILD_TYPE MATCHES Debug)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none -g -gdwarf -Og -fbacktrace -fcheck=all")
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ")
      SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ")
    else ()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  ")
      SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ")
    endif ()
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "PGI")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fPIC")
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (CMAKE_BUILD_TYPE MATCHES Debug)
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -g -O0 -traceback -check all")
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -mkl")
      SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -mkl")
    else ()
      set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -O2")
      SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}  -mkl")
      SET(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -mkl")
    endif ()
  endif ()
endif()

####################################################################################
# 4. MPI
####################################################################################
if(_is_sw_platform)
  find_package(MPI)
  include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})
else()
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    find_package(MPI)
    include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    find_package(MPI)
    include_directories(SYSTEM ${MPI_Fortran_INCLUDE_DIRS})
  endif ()
endif()

####################################################################################
# 5. NETCDF
####################################################################################
if (DEFINED ENV{NETCDF} AND (NOT DEFINED ENV{NETCDF_ROOT}))
  set(ENV{NETCDF_ROOT} $ENV{NETCDF})
endif ()

if(_is_sw_platform)
  if (DEFINED ENV{NETCDF_ROOT})
    include_directories("$ENV{NETCDF_ROOT}/include")
    link_directories("$ENV{NETCDF_ROOT}/lib")
    include_directories("/usr/sw/yyzlib/hdf5-1.8.9/include")
    link_directories("/usr/sw/yyzlib/hdf5-1.8.9/lib")
  else ()
    message(FATAL_ERROR "Pkg-config could not find sw netcdf library!")
    message(FATAL_ERROR "Please set NETCDF_ROOT path or NETCDF path!")
  endif()
else()
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
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
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
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
  endif ()
endif()

if (DEFINED ENV{HDF5})
  include_directories("$ENV{HDF5}/include")
  link_directories("$ENV{HDF5}/lib")
  set(EXTERNAL_LIBS ${EXTERNAL_LIBS} hdf5 hdf5_hl z)
  message(STATUS "Add HDF5 $ENV{HDF5} into compilation process.")
endif ()

####################################################################################
# 6. MKL
####################################################################################
if(_is_sw_platform)
  if (DEFINED ENV{XMATH_PATH})

    if(_is_sw_master)
      include_directories("$ENV{XMATH_PATH}/include")
      link_directories("/usr/sw/yyzlib/lapack-3.8.0")
      link_directories("$ENV{XMATH_PATH}/")
    else()
      #include_directories("$ENV{XMATH_PATH}/include")
      #link_directories("$ENV{XMATH_PATH}/")
	  link_directories("/home/export/online1/mdt00/shisuan/swthhpca/DIDA-SHARE/3rdParty/xMath3.0")
	  message("found XMATH!")
    endif()
  else()
    message(FATAL_ERROR "could not find xmath library!")
    message(FATAL_ERROR "please set XMATH_ROOT path!")
  endif()
else()
  if (CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(MKL_LIBRARIES "lapack64")
  elseif (CMAKE_Fortran_COMPILER_ID STREQUAL "Intel")
    if (DEFINED ENV{MKLROOT})
      include_directories("$ENV{MKLROOT}/include")
      link_directories("$ENV{MKLROOT}/lib/intel64")
    elseif (DEFINED ENV{MKL_INC} AND DEFINED ENV{MKL_LIB})
      include_directories("$ENV{MKL_INC}")
      link_directories("$ENV{MKL_LIB}")
    else ()
      find_package(PkgConfig)
      if (PKG_CONFIG_FOUND)
        pkg_search_module(MKL REQUIRED mkl)
        if (MKL_FOUND)
          include_directories(${MKL_INCLUDEDIR})
          link_directories(${MKL_LIBDIR})
        else ()
          message(FATAL_ERROR "Pkg-config could not find mkl library!")
        endif ()
      else ()
        message(FATAL_ERROR "Unable to find pkg-config library!")
      endif ()
    endif ()
  endif ()
endif()

####################################################################################
# 7. include
####################################################################################
include_directories(build/)
include_directories(build/LIB/MCT)
include_directories(build/LIB/REDIS)
include_directories(build/LIB/container)
include_directories(build/LIB/fortran_datetime)
include_directories(build/LIB/fiona)
include_directories(build/LIB/flogger)
include_directories(build/LIB/string)
include_directories(build/DA/GRNG)
include_directories(build/COUPLER)
include_directories(build/DA/utils)
include_directories(build/DA/LETKF)
include_directories(build/ATM)
include_directories(build/OBS)
include_directories(LIB/mpal/include)
include_directories(DA/CNN/c_model/include)

####################################################################################
# 8. Subdirectories
####################################################################################
if(_is_sw_platform)
    #include_directories(LIB/gptl_sw/include)
    #link_directories(LIB/gptl_sw/lib)
    include_directories(LIB/gptl-5.3/include)
    link_directories(LIB/gptl-5.3/lib)
    #include_directories(LIB/gptl/include)
    #link_directories(LIB/gptl/lib)
else()
    include_directories(LIB/gptl/include)
    link_directories(LIB/gptl/lib)
    #include_directories("DA/rpe/modules")
    #link_directories("DA/rpe/lib")
endif()
link_directories(LIB/mpal/lib)
link_directories(DA/LETKF/src/letkf_half)
link_directories(/usr/sw/yyzlib/hgemm)
link_directories(DA/CNN/c_model/libs)

add_subdirectory(LIB/MCT)
add_subdirectory(LIB/REDIS)
add_subdirectory(LIB/flogger/)
add_subdirectory(LIB/container)
add_subdirectory(LIB/string)
add_subdirectory(LIB/flogging/)
add_subdirectory(LIB/fortran-datetime/)
add_subdirectory(LIB/fiona/)
add_subdirectory(LIB/unit-test/)

add_subdirectory(COUPLER)
add_subdirectory(ATM)
add_subdirectory(DA/utils)

####################################################################################
# 9. Target_link
####################################################################################
add_executable(dida.exe DRIVE/dida_drive.F90)
if(_is_sw_platform)
    if(_is_sw_master)
        target_link_libraries(dida.exe da atm coupler redis mct flogger fortran_datetime fiona fortran_string rt m z gfortran gptl unwind unwind-ptrace unwind-setjmp unwind-generic )
    else()
        target_link_libraries(dida.exe da atm coupler redis mct flogger fortran_datetime fiona fortran_string rt m_slave z gfortran gptl unwind unwind-ptrace unwind-setjmp unwind-generic)
    endif()
else()
  target_link_libraries(dida.exe da atm coupler gptl redis mct flogger fortran_datetime fiona fortran_string librpe.a)
endif()


####################################################################################
# 10. UNITTEST
####################################################################################
if(ATM_TEST)
  add_executable(gmcore_swm_driver.exe DRIVE/gmcore_swm_driver.F90)
  target_link_libraries(gmcore_swm_driver.exe atm )
  add_executable(gmcore_driver.exe DRIVE/gmcore_driver.F90)
  target_link_libraries(gmcore_driver.exe atm )
endif()

if (fiona_TEST)
  add_executable(fiona_test.exe LIB/fiona/src/fiona_test.F90)
  target_link_libraries(fiona_test.exe fiona fortran_unit_test  netcdff netcdf hdf5_hl hdf5 z )
endif ()
