cmake_minimum_required(VERSION 3.0)

SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -std=c99 -O3 -DNOUNDERSCORE")

set(CMAKE_Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}")
include_directories(hiredis-vip/)

set (source_files
     hiredis-vip/adlist.c
     hiredis-vip/async.c
     hiredis-vip/command.c
     hiredis-vip/crc16.c
     hiredis-vip/dict.c
     hiredis-vip/hiarray.c
     hiredis-vip/hircluster.c
     hiredis-vip/hiredis.c
     hiredis-vip/hiutil.c
     hiredis-vip/net.c
     hiredis-vip/read.c
     hiredis-vip/sds.c
     hiredis-vip/logger.c
     src/redis_wrapper.c
     src/sb.c
     src/redis_mod.F90
    )

add_library(redis ${source_files})
target_link_libraries(redis fortran_string)
