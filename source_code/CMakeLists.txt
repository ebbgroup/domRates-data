cmake_minimum_required(VERSION 2.6)

SET(MAJOR_VERSION 1)
SET(MINOR_VERSION 0)
SET(PATCH_VERSION 0)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
	message(STATUS "Setting build type to 'Release' as none was specified.")
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
endif()

SET(CMAKE_CXX_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-default-inline"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
SET(CMAKE_CXX_FLAGS_PROFILE
    "-O3 -pg"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
SET(CMAKE_C_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-default-inline"
    CACHE STRING "Flags used by the C compiler during coverage builds."
    FORCE )
SET(CMAKE_EXE_LINKER_FLAGS_COVERAGE
    "--coverage"
    CACHE STRING "Flags used for linking binaries during coverage builds."
    FORCE )
SET(CMAKE_SHARED_LINKER_FLAGS_COVERAGE
    "--coverage"
    CACHE STRING "Flags used by the shared libraries linker during coverage builds."
    FORCE )
MARK_AS_ADVANCED(
    CMAKE_CXX_FLAGS_COVERAGE
    CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )


project (domRates C CXX)

if (UNIX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunreachable-code -Wall -pedantic -Wextra -Wunused -lpthread ")
	set(Boost_USE_STATIC_LIBS OFF)
	set(Boost_USE_MULTITHREADED OFF)
	if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
		if(CMAKE_CXX_COMPILER_VERSION VERSION_LESS "4.7")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
		endif()
		if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.6")
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
		endif()
	endif()
endif(UNIX)



################################################
#         find necessary libraries             #
################################################

# boost
if (WITH_UNIT_TEST)
	FIND_PACKAGE(Boost 1.54 COMPONENTS system program_options iostreams filesystem unit_test_framework REQUIRED)
else (WITH_UNIT_TEST)
	FIND_PACKAGE(Boost 1.54 COMPONENTS system program_options iostreams filesystem REQUIRED)
endif (WITH_UNIT_TEST)

# boost
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
ADD_DEFINITIONS( "-DHAS_BOOST" )

# OpenMP
cmake_policy(SET CMP0012 NEW)
find_package(OpenMP 4.0)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()



################################################
#                   domRates                   #
################################################


include_directories(${CMAKE_CURRENT_BINARY_DIR})

set (BSDL_PATH ./libs/BioSeqDataLib/src/) 
set (BSDL_src ${BSDL_PATH}/external/Input.cpp ${BSDL_PATH}/external/Output.cpp ${BSDL_PATH}/domain/Domain.cpp ${BSDL_PATH}/domain/PfamDomain.cpp ${BSDL_PATH}/domain/DomainExt.cpp ${BSDL_PATH}/domain/SFDomain.cpp ${BSDL_PATH}/domain/DomainArrangement.cpp ${BSDL_PATH}/domain/DomainArrangementSet.cpp ${BSDL_PATH}/phylogeny/PhylogeneticTree.cpp ${BSDL_PATH}/phylogeny/fitch.cpp ${BSDL_PATH}/phylogeny/dollo.cpp ${BSDL_PATH}/utility/stringHelpers.cpp)
set (DOMRA_PATH ./src/)
set (DOMRA_src ${DOMRA_PATH}/domRates.cpp)

SET(domRates_src src/domRates_main.cpp ${BSDL_src} ${DOMRA_src})
SET(domRates_exe domRates)
ADD_EXECUTABLE(${domRates_exe} ${domRates_src})
target_link_libraries(${domRates_exe} 
	${Boost_LIBRARIES}
)

configure_file (
		"${PROJECT_SOURCE_DIR}/src/version.hpp"
		"${PROJECT_SOURCE_DIR}/src/version.hpp"
)

if (WITH_UNIT_TEST)
	add_subdirectory (tests)
	enable_testing ()
	add_test (NAME UnitTest COMMAND unit_tests)
endif ()
