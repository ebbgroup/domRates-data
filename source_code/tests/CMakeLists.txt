
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED OFF)

FUNCTION(PREPEND var prefix)
    SET(listVar "")
    FOREACH(f ${ARGN})
        LIST(APPEND listVar "${prefix}${f}")
    ENDFOREACH(f)
    SET(${var} "${listVar}" PARENT_SCOPE)
ENDFUNCTION(PREPEND)

SET(tests_src ./unitTests/unit_tests.cpp ../src/domRates.cpp ../libs/BioSeqDataLib/src/phylogeny/PhylogeneticTree.cpp ../libs/BioSeqDataLib/src/domain/Domain.cpp ../libs/BioSeqDataLib/src/domain/DomainArrangement.cpp ../libs/BioSeqDataLib/src/domain/DomainArrangementSet.cpp ../libs/BioSeqDataLib/src/domain/DomainExt.cpp ../libs/BioSeqDataLib/src/domain/PfamDomain.cpp ../libs/BioSeqDataLib/src/domain/SFDomain.cpp ../libs/BioSeqDataLib/src/phylogeny/fitch.cpp ../libs/BioSeqDataLib/src/phylogeny/dollo.cpp ../libs/BioSeqDataLib/src/external/Input.cpp ../libs/BioSeqDataLib/src/external/Output.cpp ../libs/BioSeqDataLib/src/utility/stringHelpers.cpp)
SET(tests_exe unit_tests)
ADD_EXECUTABLE(${tests_exe} ${tests_src})
target_link_libraries(${tests_exe}
        ${Boost_LIBRARIES}
        )
