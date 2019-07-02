
#include "common.h"
#include "catch.hpp"
#include "../Parser.h"
#include <stdio.h>
#include <string>
#include <iostream>

// system(diff) in this testing case is system dependent
// this is ok for testing purposes

int returnDiff(string expected, string result) {
    expected = PATH + "/test_files/" + expected;
    string command = "diff " + expected + " " + result + " >/dev/null";
    return system(command.c_str());
}

TEST_CASE("Test whole program with different settings") {
    Parser fileParser;
    fileParser.setWindowLegth(10);

    string inputFile = PATH + "/test_files/test_1.ali";
    string output = PATH + "/test_files/test_result";
    ScoreMatrix * scoreMatrix = new ScoreMatrix();
    scoreMatrix->loadFromFile(PATH + "/test_files/blosum62_1.csv");
    Kernel * triangularkernel = new TriangularKernel();
    fileParser.setScoringMatrix(scoreMatrix);
    fileParser.setKernel(triangularkernel);

    SECTION("All") {
        freopen(inputFile.c_str(), "r", stdin);
        fileParser.setMinExonScore(-999999);
        fileParser.parse(output);
        int result = returnDiff("test_parser.gff", output);
        CHECK(result == 0);
    }

    std::cin.clear();

    SECTION("Exon filter 25") {
        freopen(inputFile.c_str(), "r", stdin);
        fileParser.setMinExonScore(25);
        fileParser.parse(output);
        int result = returnDiff("test_parser_eScore_25.gff", output);
        CHECK(result == 0);
    }

    delete scoreMatrix;
    delete triangularkernel;
    remove(output.c_str());
}
