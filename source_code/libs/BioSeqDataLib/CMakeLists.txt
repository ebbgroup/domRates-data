cmake_minimum_required(VERSION 2.8)

# onec cmake 3 will be used
#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-std=c++11)

project (BioSeqDataLib_TESTS CXX)
SET(MAJOR_VERSION 0)
SET(MINOR_VERSION 9)
SET(PATCH_VERSION 9)

SET(CMAKE_CXX_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-default-inline -pedantic -Wall -Wextra" # -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Wno-unused"
    CACHE STRING "Flags used by the C++ compiler during coverage builds."
    FORCE )
SET(CMAKE_C_FLAGS_COVERAGE
    "-g -O0 --coverage -fprofile-arcs -ftest-coverage -fno-inline -fno-inline-small-functions -fno-default-inline -pedantic -Wall -Wextra" # -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wshadow -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Wno-unused"
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
  #  CMAKE_C_FLAGS_COVERAGE
    CMAKE_EXE_LINKER_FLAGS_COVERAGE
    CMAKE_SHARED_LINKER_FLAGS_COVERAGE )



if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weverything -Wno-padded -Wno-c++98-compat ")
endif ()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -Wextra -Wunused")# -Wcast-align -Wcast-qual -Wctor-dtor-privacy -Wdisabled-optimization -Wformat=2 -Winit-self -Wlogical-op -Wmissing-declarations -Wmissing-include-dirs -Wnoexcept -Wold-style-cast -Woverloaded-virtual -Wredundant-decls -Wsign-conversion -Wsign-promo -Wstrict-null-sentinel -Wstrict-overflow=5 -Wswitch-default -Wundef -Werror -Wno-unused")
endif ()
#message(${CMAKE_CXX_FLAGS})
#-Wshadow


if( NOT CMAKE_BUILD_TYPE )
  set( CMAKE_BUILD_TYPE Release CACHE STRING
       "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel." FORCE )
endif()



FUNCTION(PREPEND var prefix)
   SET(listVar "")
   FOREACH(f ${ARGN})
      LIST(APPEND listVar "${prefix}/${f}")
   ENDFOREACH(f)
   SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)



# To run in debug mode -DCMAKE_BUILD_TYPE=Debug

if (WITH_UNIT_TEST)
	FIND_PACKAGE(Boost 1.49 COMPONENTS system filesystem iostreams unit_test_framework REQUIRED)
else (WITH_UNIT_TEST)
	FIND_PACKAGE(Boost 1.49 COMPONENTS system filesystem iostreams REQUIRED)
endif(WITH_UNIT_TEST)
INCLUDE_DIRECTORIES(SYSTEM ${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
ADD_DEFINITIONS( "-DHAS_BOOST" )


# The annotation module
set(annotationCPP Feature.cpp BlastHit.cpp FeatureSet.cpp OrthologySet.cpp)
PREPEND(annotationCPP "${CMAKE_CURRENT_SOURCE_DIR}/src/annotation" ${annotationCPP})

# The domain module
set(domainCPP Domain.cpp DomainExt.cpp PfamDomain.cpp SFDomain.cpp DomainArrangementSet.cpp DomainArrangement.cpp )
PREPEND(domainCPP "${CMAKE_CURRENT_SOURCE_DIR}/src/domain" ${domainCPP})

# external_interfaces module
set(external_interfacesCPP runPhylo.cpp runDBsearch.cpp domainProgs.cpp)
PREPEND(external_interfacesCPP "${CMAKE_CURRENT_SOURCE_DIR}/src/external_interfaces" ${external_interfacesCPP})

# phylogeny module
set(phylogenyCPP PhylogeneticTree.cpp fitch.cpp dollo.cpp)
PREPEND(phylogenyCPP "${CMAKE_CURRENT_SOURCE_DIR}/src/phylogeny" ${phylogenyCPP})

# utility module
set(utilityCPP TwoValues.cpp DSM.cpp stringHelpers.cpp properties.cpp utility.cpp Settings.cpp)
PREPEND(utilityCPP "${CMAKE_CURRENT_SOURCE_DIR}/src/utility" ${utilityCPP})

#external module
set(externalCPP Input.cpp Output.cpp)
PREPEND(externalCPP "${CMAKE_CURRENT_SOURCE_DIR}/src/external" ${externalCPP})


set(SOURCE_FILES ${annotationCPP} ${external_interfacesCPP} ${externalCPP} ${phylogenyCPP} ${utilityCPP} ${domainCPP})
IF (CURL_FOUND)
	set(SOURCE_FILES ${SOURCE_FILES} ${webserviceCPP})
ENDIF(CURL_FOUND)


add_library(BioSeqDataLib SHARED ${SOURCE_FILES})
set_target_properties (BioSeqDataLib PROPERTIES VERSION ${MAJOR_VERSION}.${MINOR_VERSION}.${PATCH_VERSION})
target_link_libraries(BioSeqDataLib ${Boost_LIBRARIES}  ${CURL_LIBRARIES})


INSTALL(TARGETS BioSeqDataLib
	DESTINATION "lib"
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/src/ DESTINATION include/BioSeqDataLib
	FILES_MATCHING PATTERN "*.hpp"
	PATTERN ${CMAKE_CURRENT_SOURCE_DIR}/src/documentation.hpp EXCLUDE
)

if (WITH_UNIT_TEST)
	add_subdirectory (tests)
	enable_testing ()
	add_test (NAME domainTest COMMAND domain_tests)
	add_test (NAME SequenceTest COMMAND sequence_tests)
	add_test (NAME PhylogenyTest COMMAND phylogeny_tests)
	add_test (NAME AnnotationTest COMMAND annotation_tests)
	add_test (NAME AlignTest COMMAND align_tests)
	add_test (NAME external_interfacesTest COMMAND external_interfaces_tests)
	add_test (NAME UtilityTest COMMAND utility_tests)
endif ()
