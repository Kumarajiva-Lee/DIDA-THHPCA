cmake_minimum_required(VERSION 3.0)

include_directories($PROJECT_SOURCE_DIR/build/LIB/string)

set (source_files
     src/module_coupler.F90
     src/module_coupler_config.F90
     src/module_coupler_log.F90
    )

add_library(coupler ${source_files})
target_link_libraries(coupler redis mct fortran_string)
