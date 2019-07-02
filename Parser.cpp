#include "Parser.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

Parser::Parser() {
    windowLength = 10;
    scoreMatrix = NULL;
}

int Parser::parse(string outputFile) {
    // Restart output file
    ofstream ofs(outputFile.c_str());
    ofs.close();

    int status = parseNext();

    while (status != NO_MORE_ALIGNMENTS) {
        alignment.scoreHints(windowLength, scoreMatrix, kernel);
        alignment.printHints(outputFile, minExonScore);
        status = parseNext();
    }

    return READ_SUCCESS;
}

int Parser::parseNext() {
    string line;
    while (getline(cin, line)) {
        if (line.substr(0,1) == ">") {
            return alignment.parse(cin, line);
        }
    }

    return NO_MORE_ALIGNMENTS;
}

double Parser::maxScore() {
    return scoreMatrix->getMaxScore();
}

void Parser::setScoringCombination(int combination) {
    scoreCombination = combination;
}

void Parser::setWindowLegth(int length) {
    windowLength = length;
}

void Parser::setScoringMatrix(const ScoreMatrix * scoreMatrix) {
    this->scoreMatrix = scoreMatrix;
}

void Parser::setKernel(Kernel* kernel) {
    this->kernel = kernel;
}

void Parser::setMinExonScore(double minExonScore) {
    this->minExonScore = minExonScore;
}
