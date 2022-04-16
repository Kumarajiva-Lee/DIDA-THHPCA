cmake_minimum_required(VERSION 3.0)

set (source_files
     src/module_da.F90
     src/module_da_config.F90
     src/module_da_coupler.F90
     src/module_da_parallel.F90
    )

add_library(endnn ${source_files})
target_link_libraries(endnn coupler redis mct)
