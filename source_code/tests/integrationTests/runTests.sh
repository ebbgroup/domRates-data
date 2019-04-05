#!/usr/bin/env bats

@test "DomRates IntegrationTest" {

	# test run without detailed output (-d parameter)
	run ../../build/domRates -t ../data/test_tree.nwk -a ../data/ -g OG -o results/test_out.txt -s results/test_out_stats.txt
	[ $status == 0 ]
	run diff -I '# DomRates' results/compare_out.txt results/test_out.txt
	[ $status == 0 ]
	run diff results/compare_out_stats.txt results/test_out_stats.txt
	[ $status == 0 ]
	run diff results/compare_out_stats_epd.txt results/test_out_stats_epd.txt
	[ $status == 0 ]

    	# run parallelization test with 4 threads
    	#run ../../build/domRates -t ../data/test_tree.nwk -a ../data -g OG -o results/test_out.txt -p 4
	#[ $status == 0 ]
	#run diff <(results/compare_output.txt) <(results/test_out.txt)
	#[ $status == 0 ]

	# test run with detailed output (-d parameter)
	run ../../build/domRates -t ../data/test_tree.nwk -a ../data/ -g OG -o results/test_ident_out.txt -s results/test_ident_out_stats.txt -d
        [ $status == 0 ]
        run diff -I '# DomRates' results/compare_ident_out.txt results/test_ident_out.txt
	[ $status == 0 ]
        run diff results/compare_ident_out_stats.txt results/test_ident_out_stats.txt
	[ $status == 0 ]
	run diff results/compare_ident_out_stats_epd.txt results/test_ident_out_stats_epd.txt
        [ $status == 0 ]

	# test run with detailed output (-d parameter) and statistics for LCA node of A and B (-n parameter)
	run ../../build/domRates -t ../data/test_tree.nwk -a ../data/ -g OG -o results/test_nident_out.txt -s results/test_nident_out_stats.txt -d -n A:B
	[ $status == 0 ]
	run diff -I '# DomRates' results/compare_nident_out.txt results/test_nident_out.txt
	[ $status == 0 ]
	run diff results/compare_nident_out_stats.txt results/test_nident_out_stats.txt
	[ $status == 0 ]
	run diff results/compare_nident_out_stats_epd.txt results/test_nident_out_stats_epd.txt
	[ $status == 0 ]

	
	# test run without output file (rates printed in console)
	run ../../build/domRates -t ../data/test_tree.nwk -a ../data/ -g OG -d
	[ $status == 0 ]
	run diff <(echo $output) <(cat results/test_ident_out.txt)

	run ../../build/domRates -h
	[ $status == 0 ]

	rm results/test_out.txt
	rm results/test_out_stats.txt
	rm results/test_out_stats_epd.txt

	rm results/test_ident_out.txt
        rm results/test_ident_out_stats.txt
        rm results/test_ident_out_stats_epd.txt

	rm results/test_nident_out.txt
	rm results/test_nident_out_stats.txt
	rm results/test_nident_out_stats_epd.txt
}

