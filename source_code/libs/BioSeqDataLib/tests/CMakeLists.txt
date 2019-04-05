
SET(sequence_tests_src ./sequence/sequence_tests.cpp)
SET(sequence_tests_exe sequence_tests)
ADD_EXECUTABLE(${sequence_tests_exe} ${sequence_tests_src})
target_link_libraries(${sequence_tests_exe}
	${Boost_LIBRARIES} BioSeqDataLib
)

SET(annotation_tests_src ./annotation/annotation_tests.cpp)
SET(annotation_tests_exe annotation_tests)
ADD_EXECUTABLE(${annotation_tests_exe} ${annotation_tests_src})
target_link_libraries(${annotation_tests_exe}
	${Boost_LIBRARIES} BioSeqDataLib
)

SET(align_tests_src ./align/alignTests.cpp)
SET(align_tests_exe align_tests)
ADD_EXECUTABLE(${align_tests_exe} ${align_tests_src})
target_link_libraries(${align_tests_exe}
	${Boost_LIBRARIES} BioSeqDataLib
)

SET(domain_tests_src ./domain/domain_tests.cpp)
SET(domain_tests_exe domain_tests)
ADD_EXECUTABLE(${domain_tests_exe} ${domain_tests_src})
target_link_libraries(${domain_tests_exe}
	${Boost_LIBRARIES} BioSeqDataLib
)

SET(utility_tests_src ./utility/utility_tests.cpp)
SET(utility_tests_exe utility_tests)
ADD_EXECUTABLE(${utility_tests_exe} ${utility_tests_src})
target_link_libraries(${utility_tests_exe}
	${Boost_LIBRARIES} BioSeqDataLib
)

SET(phylogeny_tests_src ./phylogeny/phylogeny_tests.cpp)
SET(phylogeny_tests_exe phylogeny_tests)
ADD_EXECUTABLE(${phylogeny_tests_exe} ${phylogeny_tests_src})
target_link_libraries(${phylogeny_tests_exe}
	${Boost_LIBRARIES} BioSeqDataLib
)

SET(external_interfaces_tests_src ./external_interfaces/external_interfaces_tests.cpp)
SET(external_interfaces_tests_exe external_interfaces_tests)
ADD_EXECUTABLE(${external_interfaces_tests_exe} ${external_interfaces_tests_src})
target_link_libraries(${external_interfaces_tests_exe}
	${Boost_LIBRARIES} BioSeqDataLib
)
