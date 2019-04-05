
#ifndef DOMRATES_FILEHANDLING_TEST_HPP
#define DOMRATES_FILEHANDLING_TEST_HPP

#include <boost/test/unit_test.hpp>
#include <iostream>
#include <string>

#include <boost/filesystem.hpp>

#include "../../src/domRates.hpp"

BOOST_AUTO_TEST_SUITE(Filehandling_Test)

    BOOST_AUTO_TEST_CASE(alterFilename)
    {
        std::string appdx = "_epd";

        boost::filesystem::path file1 = alter_filename("file_1", appdx);
        BOOST_CHECK_EQUAL(file1.string(), "file_1_epd.txt");

        fs::path file2 = alter_filename("/path/to/smth/file_2", appdx);
        BOOST_CHECK_EQUAL(file2.string(), "/path/to/smth/file_2_epd.txt");

        fs::path file3 = alter_filename("file_3.txt", appdx);
        BOOST_CHECK_EQUAL(file3.string(), "file_3_epd.txt");

        fs::path file4 = alter_filename("/path/to/smth/file_4.txt", appdx);
        BOOST_CHECK_EQUAL(file4.string(), "/path/to/smth/file_4_epd.txt");

        fs::path file5 = alter_filename("file_5.out", appdx);
        BOOST_CHECK_EQUAL(file5.string(), "file_5_epd.out");

        fs::path file6 = alter_filename("/path/to/smth/file_6.out", appdx);
        BOOST_CHECK_EQUAL(file6.string(), "/path/to/smth/file_6_epd.out");

        fs::path file7 = alter_filename("file.with.many.ext.out", appdx);
        BOOST_CHECK_EQUAL(file7.string(), "file.with.many.ext_epd.out");

        fs::path file8 = alter_filename("/path/to/smth/file.with.many.ext.out", appdx);
        BOOST_CHECK_EQUAL(file8.string(), "/path/to/smth/file.with.many.ext_epd.out");
    }

BOOST_AUTO_TEST_SUITE_END()

#endif //DOMRATES_FILEHANDLING_TEST_HPP
