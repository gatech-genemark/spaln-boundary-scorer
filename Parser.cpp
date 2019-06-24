#include "Parser.h"
#include <string>
#include <fstream>
#include <iostream>

using namespace std;

//const char Parser::PAIR_SEPARATOR[] = "************************************"
//        "************************************";

Parser::Parser() {
    scoreCombination = BOUNDARIES_SUMMED;
    windowLength = 20;
    printAll = false;
    scoreMatrix = NULL;
}

int Parser::parse(string intputFile, string outputFile) {
    ofstream ofs(outputFile.c_str());
    ofs.close();
    inputStream.open(intputFile.c_str());
    if (!inputStream) {
        cerr << "error: Failed to open file \"" << intputFile << "\"" << endl;
        return OPEN_FAIL;
    }
    IntronStorage storage;
    int status = parseNext();

    while (status != NO_MORE_ALIGNMENTS) {
        //cout << alignment.getGene() << " " << alignment.getProtein() << endl;
        if (scoreCombination == BOUNDARIES_SUMMED) {
            alignment.scoreIntrons(windowLength, false, scoreMatrix, kernel);
        } else {
            alignment.scoreIntrons(windowLength, true, scoreMatrix, kernel);
        }
        alignment.storeIntrons(outputFile, storage);
        status = parseNext();
    }
    storage.normalizeScores(0, maxScore());
    storage.printIntrons(outputFile, printAll);
    inputStream.close();
    return READ_SUCCESS;
}

int Parser::parseNext() {
    string line;
    while (getline(inputStream, line)) {
        if (line.substr(0,1) == ">") {
            return alignment.parse(inputStream, line);
        }
    }

    return NO_MORE_ALIGNMENTS;
}

double Parser::maxScore() {
    if (scoreMatrix == NULL) {
        return COMPLETE_MATCH;
    } else {
        return scoreMatrix->getMaxScore();
    }
}

void Parser::setScoringCombination(int combination) {
    scoreCombination = combination;
}

void Parser::setWindowLegth(int length) {
    windowLength = length;
}

void Parser::printAllSites() {
    printAll = true;
}

void Parser::printCanonicalsites() {
    printAll = false;
}

void Parser::setScoringMatrix(const ScoreMatrix * scoreMatrix) {
    this->scoreMatrix = scoreMatrix;
}

void Parser::setKernel(Kernel* kernel) {
    this->kernel = kernel;
}
