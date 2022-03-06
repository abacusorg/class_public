cmake_minimum_required(VERSION 3.11)
project(class)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
# set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fsanitize=address")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

# Options
option(ENABLE_TESTS "compile test executables" OFF)
option(ENABLE_HYREC "compile hyrec in ./external/HyRec2020 directory" ON)
option(ENABLE_RECFAST "compile recfast in ./external/RecfastCLASS directory" ON)
option(ENABLE_HEATING "compile heating in ./external/heating directory" ON)

# OpenMP
find_package(OpenMP)

# Sources
set(TOOLS_SRC 
    tools/growTable.c 
    tools/dei_rkck.c 
    tools/sparse.c 
    tools/evolver_rkck.c  
    tools/evolver_ndf15.c 
    tools/arrays.c 
    tools/parser.c 
    tools/quadrature.c 
    tools/hyperspherical.c 
    tools/common.c 
    tools/trigonometric_integrals.c
)

set(SOURCE_SRC
    source/input.c 
    source/output.c
    source/background.c 
    source/thermodynamics.c 
    source/perturbations.c 
    source/primordial.c 
    source/fourier.c
    source/transfer.c 
    source/harmonic.c
    source/lensing.c
    source/distortions.c
)

set(CPP_SRC 
    cpp/ClassEngine.hh cpp/ClassEngine.cc 
    cpp/Engine.hh cpp/Engine.cc
)

set(TESTS
    test_loops
    test_loops_omp
    test_spectra
    test_transfer
    test_nonlinear
    test_perturbations
    test_thermodynamics
    test_background
    test_hyperspherical
)


set(EXTERNAL_SRC )
set(EXTERNAL_INCLUDE )

if(ENABLE_HYREC)
    list(APPEND EXTERNAL_SRC
         external/HyRec2020/hyrectools.c
         external/HyRec2020/helium.c
         external/HyRec2020/hydrogen.c
         external/HyRec2020/history.c
         external/HyRec2020/wrap_hyrec.c
         external/HyRec2020/energy_injection.c
    )
    list(APPEND EXTERNAL_INCLUDE external/HyRec2020)
endif()

if(ENABLE_RECFAST)
    list(APPEND EXTERNAL_SRC
        external/RecfastCLASS/wrap_recfast.c
    )
    list(APPEND EXTERNAL_INCLUDE external/RecfastCLASS)
endif()

if(ENABLE_HEATING)
    list(APPEND EXTERNAL_SRC
        external/heating/injection.c
        external/heating/noninjection.c
    )
    list(APPEND EXTERNAL_INCLUDE external/heating)
endif()

## Main library
add_library(libclass ${TOOLS_SRC} ${SOURCE_SRC} ${EXTERNAL_SRC})
add_library(class::libclass ALIAS libclass)
set_target_properties(libclass PROPERTIES OUTPUT_NAME class)
set_target_properties(libclass PROPERTIES C_STANDARD 11)
set_property(TARGET libclass PROPERTY POSITION_INDEPENDENT_CODE ON)
# target_compile_options(libclass PRIVATE "-O4" "-ffast-math" "-g")
target_compile_definitions(libclass PUBLIC "__CLASSDIR__=\"${CMAKE_CURRENT_SOURCE_DIR}\"")
target_include_directories(libclass PUBLIC include)
target_include_directories(libclass PUBLIC ${EXTERNAL_INCLUDE})
target_link_libraries(libclass PRIVATE m)
if(OpenMP_FOUND)
    target_link_libraries(libclass PRIVATE OpenMP::OpenMP_C)
endif()
if(ENABLE_HYREC)
    target_compile_definitions(libclass PUBLIC "HYREC")
    target_include_directories(libclass PUBLIC ${CMAKE_CURRENT_LIST_DIR}/external/HyRec2020)
endif()
if(ENABLE_RECFAST)
    target_compile_definitions(libclass PUBLIC "RECFAST")
    target_include_directories(libclass PUBLIC ${CMAKE_CURRENT_LIST_DIR}/external/RecfastCLASS)
endif()
if(ENABLE_HEATING)
    target_compile_definitions(libclass PUBLIC "HEATING")
    target_include_directories(libclass PUBLIC ${CMAKE_CURRENT_LIST_DIR}/external/heating)
endif()

## Main Executable
add_executable(class main/class.c)
target_link_libraries(class PRIVATE libclass)

## Tests
if(ENABLE_TESTS)
    foreach(test ${TESTS})
        add_executable(${test} "test/${test}.c")
        # target_compile_options(${test} PRIVATE "-O4" "-ffast-math" "-g")
        target_link_libraries(${test} PRIVATE libclass)
        if(OpenMP_FOUND)
            target_link_libraries(${test} PRIVATE OpenMP::OpenMP_C)
        endif()
    endforeach()
endif()

## C++ library
add_library(libclass_cpp ${CPP_SRC})
add_library(class::libclass_cpp ALIAS libclass_cpp)
set_target_properties(libclass_cpp PROPERTIES CXX_STANDARD 14)
target_link_libraries(libclass_cpp PUBLIC libclass)
target_include_directories(libclass_cpp PUBLIC cpp)
set_target_properties(libclass_cpp PROPERTIES OUTPUT_NAME class_cpp)



## Python (classy)
## TODO ##
