cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(PETSIRDUseCaseTemplate VERSION 0.7.2)

set(CMAKE_CXX_STANDARD 17)

#Set the build type to Release if not specified
if (NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING
            "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
            FORCE)
endif ()

if(WIN32)
  add_compile_options(/W3 /WX)
else()
  add_compile_options(-Wall -Wextra -pedantic)
endif()

set(PETSIRD_dir ../PETSIRD/cpp/)
add_subdirectory(${PETSIRD_dir} PETSIRD)


add_executable(petsird_castor_converter petsird_castor_converter.cpp)
# target_include_directories(petsird_castor_converter PRIVATE "${PETSIRD_dir}/")
# needed for generated PETSIRD files
# target_include_directories(petsird_castor_converter PRIVATE "${PETSIRD_dir}/../")
# target_include_directories(petsird_castor_converter PRIVATE "${PETSIRD_dir}/../helpers/include")
# target_link_libraries(petsird_castor_converter petsird_generated)
target_link_libraries(petsird_castor_converter PUBLIC petsird_helpers)

# # install(TARGETS my_prog DESTINATION bin)
