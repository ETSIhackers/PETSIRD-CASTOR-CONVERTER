cmake_minimum_required(VERSION 3.12.0) # older would work, but could give warnings on policy CMP0074
project(PETSIRDUseCaseTemplate VERSION 0.7.2)

set(CMAKE_CXX_STANDARD 20)
<<<<<<< HEAD
=======
set(CMAKE_CXX_STANDARD_REQUIRED ON)

>>>>>>> change_centroid_to_ADP

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

<<<<<<< HEAD
#set(PETSIRD_dir ../PETSIRD/cpp/)
#add_subdirectory(${PETSIRD_dir} PETSIRD)
=======
>>>>>>> change_centroid_to_ADP

find_package(PETSIRD 0.13 REQUIRED CONFIG)
#target_include_directories(petsird_castor_converter PUBLIC PETSIRD::helpers)

find_package(PETSIRD 0.13 REQUIRED CONFIG)
#target_include_directories(petsird_castor_converter PUBLIC PETSIRD::helpers)

add_executable(petsird_castor_converter petsird_castor_converter.cpp)
# target_include_directories(petsird_castor_converter PRIVATE "${PETSIRD_dir}/")
# needed for generated PETSIRD files
# target_include_directories(petsird_castor_converter PRIVATE "${PETSIRD_dir}/../")
# target_include_directories(petsird_castor_converter PRIVATE "${PETSIRD_dir}/../helpers/include")
# target_link_libraries(petsird_castor_converter petsird_generated)
<<<<<<< HEAD

#target_link_libraries(petsird_castor_converter PUBLIC petsird_helpers)
=======
>>>>>>> change_centroid_to_ADP
target_link_libraries(petsird_castor_converter PUBLIC PETSIRD::helpers)

# # install(TARGETS my_prog DESTINATION bin)
