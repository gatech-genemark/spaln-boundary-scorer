#include "Alignment.h"
#include "Parser.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <cmath>
#include <ctype.h>

using namespace std;

Alignment::Alignment() {
    pairs.reserve(N);
}

void Alignment::clear() {
    index = 0;
    insideIntron = false;
    donorFlag = false;
    introns.clear();
    alignmentStart = 0;
}

int Alignment::parse(ifstream& inputStream, string headerLine) {
    clear();
    vector<string> blockLines(BLOCK_ITEMS_CNT);
    string line;

    // Read header
    int status = parseHeader(headerLine);
    if (status != READ_SUCCESS) {
        cerr << "error: Invalid alignment header " << endl;
        return status;
    }
    // Get to the alignment
    getline(inputStream, line);
    while (line.substr(0,9) != "ALIGNMENT") {
        if (!getline(inputStream, line)) {
            cerr << "error: Alignment is missing after header " << endl;
            return FORMAT_FAIL;
        }
    }
    if (!getline(inputStream, line) || line != "") {
        cerr << "error: Empty line expected after ALIGNMENT keyword" << endl;
        return FORMAT_FAIL;
    }

    // Read content
    for (i = 0; i < BLOCK_ITEMS_CNT; i++) {
        if (getline(inputStream, line) && !line.empty() && line.size() > BLOCK_OFFSET) {
            blockLines[i] = line;
            // Get alignment start number
            if (i == 1) {
                alignmentStart = atoi(line.substr(0, BLOCK_OFFSET).c_str());
                if (alignmentStart <= 0) {
                    cerr << "error: Could not read alignment start position" << endl;
                    return FORMAT_FAIL;
                }
                realPositionCounter = alignmentStart;
            }
        } else {
            printLineError();
            return FORMAT_FAIL;
        }
    }

    blockLength = blockLines[1].find(" ", BLOCK_OFFSET) - BLOCK_OFFSET;
    for (i = 0; i < BLOCK_ITEMS_CNT; i++) {
        blockLines[i] = blockLines[i].substr(BLOCK_OFFSET, blockLength);
        if (blockLines[i].size() != blockLength) {
            printLineError();
            return FORMAT_FAIL;
        }
    }

    parseBlock(blockLines);

    return READ_SUCCESS;
}

int Alignment::parseHeader(string headerLine) {
    bool geneParsed = false;
    for (unsigned int i = 0; i < headerLine.size(); i++) {
        if (headerLine[i] == '>') {
            if (!geneParsed) {
                int length = headerLine.find(" ", i + 1) - i - 1;
                gene = headerLine.substr(i + 1, length);
                geneParsed = true;
            } else {
                int length = headerLine.find(" ", i + 1) - i - 1;
                protein = headerLine.substr(i + 1, length);
                return READ_SUCCESS;
            }
        }
    }
    return FORMAT_FAIL;
}

void Alignment::parseBlock(const vector<string>& lines) {
    // Parse individual pairs
    for (unsigned i = 0; i < lines[0].size(); i++) {
        AlignedPair pair(lines[0][i], lines[1][i], lines[2][i]);

        checkForIntron(pair);

        if (pair.nucleotide != '-') {
            pair.realPosition = realPositionCounter++;
        }

        // Reuse space if possible
        if ((int) pairs.size() <= index) {
            pairs.push_back(pair);
        } else {
            pairs[index] = pair;
        }
        index++;
    }
    assignCodonPhases();
}

void Alignment::printLineError() {
    cerr << "warning: error in alignment " << gene << "-" << protein;
    cerr << ": corrupted alignment - wrong line length.";
    cerr << " The rest of this alignment is skipped." << endl;
}

bool Alignment::gapOrAA(char a) {
    if ((a >= 'A' && a <= 'Z') || a == '-') {
        return true;
    }
    return false;
}

void Alignment::assignCodonPhases() {
    for (unsigned int i = 0; i < blockLength; i++) {
        if (pairs[i].translatedCodon == ' ' && pairs[i].type == 'e') {
            if (i == 0) {
                pairs[i].translatedCodon = '1';
            } else if (i == blockLength -1) {
                pairs[i].translatedCodon = '3';
            } else if (gapOrAA(pairs[i + 1].translatedCodon)) {
                pairs[i].translatedCodon = '1';
            } else if (gapOrAA(pairs[i - 1].translatedCodon)) {
                pairs[i].translatedCodon = '3';
            } else if (pairs[i + 1].type == 'i') {
                pairs[i].translatedCodon = '1';
            } else if (pairs[i - 1].type == 'i') {
                pairs[i].translatedCodon = '3';
            }
        }

        if (pairs[i].protein == ' ' && pairs[i].type == 'e') {
            if (i == 0) {
                pairs[i].protein = '1';
            } else if (i == blockLength -1) {
                pairs[i].protein = '3';
            } else if (gapOrAA(pairs[i + 1].protein)) {
                pairs[i].protein = '1';
            } else if (gapOrAA(pairs[i - 1].protein)) {
                pairs[i].protein = '3';
            } else if (pairs[i + 1].type == 'i') {
                pairs[i].protein = '1';
            } else if (pairs[i - 1].type == 'i') {
                pairs[i].protein = '3';
            }
        }
    }
}

void Alignment::checkForIntron(AlignedPair& pair) {
    if (donorFlag) {
        introns.back().donor[1] = pair.nucleotide;
        donorFlag = false;
        // If the read is at donor position, there is nothing else to check for
        return;
    }
    if (!insideIntron && pair.type == 'i') { // intron start
        Intron i;
        i.start = index;
        i.donor[0] = pair.nucleotide;
        introns.push_back(i);
        insideIntron = true;
        donorFlag = true;
    } else if (insideIntron && pair.type != 'i') { // intron end
        introns.back().end = index - 1;
        introns.back().acceptor[0] = pairs[index - 2].nucleotide;
        introns.back().acceptor[1] = pairs[index - 1].nucleotide;
        insideIntron = false;
       introns.back().complete = true;
    }
}

void Alignment::storeIntrons(IntronStorage& storage) {
    for (unsigned int i = 0; i < introns.size(); i++) {
        if (!introns[i].complete) {
            continue;
        }
        string spliceSites(introns[i].donor, 2);
        spliceSites.append("_");
        spliceSites.append(introns[i].acceptor, 2);

        if (!introns[i].scoreSet) {
            introns[i].score = 0;
        }

        storage.storeIntron(protein, gene,
                pairs[introns[i].start].realPosition,
                pairs[introns[i].end].realPosition,
                '+', spliceSites, introns[i].score, i + 1);
    }
}

void Alignment::print(ostream& os) {
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].translatedCodon;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].nucleotide;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].protein;
    }
    os << endl;
    for (unsigned int i = 0; i < blockLength; i++) {
        os << pairs[i].type;
    }
    os << endl;
}

bool Alignment::hasIntrons() {
    return introns.size() != 0;
}

string Alignment::getGene() {
    return gene;
}

string Alignment::getProtein() {
    return protein;
}

int Alignment::getLength() {
    return index;
}

Alignment::Intron::Intron() {
    scoreSet = false;
    complete = false;
}

void Alignment::scoreIntrons(int windowWidth, bool multiply,
        const ScoreMatrix * scoreMatrix, Kernel * kernel) {
    this->scoreMatrix = scoreMatrix;
    this->kernel = kernel;
    this->kernel->setWidth(windowWidth);
    for (unsigned int i = 0; i < introns.size(); i++) {
        if (introns[i].complete && !introns[i].scoreSet) {
           scoreIntron(introns[i], windowWidth, multiply);
        }
    }
}

double Alignment::scoreIntron(Intron& intron, int windowWidth, bool multiply) {
    intron.leftScore = intron.rightScore = 0;
    intron.leftWeightSum = intron.rightWeightSum = 0;
    int left, right;

    // Determine if codon is split and how
    if (pairs[intron.start - 1].protein == '3' ||
            pairs[intron.start - 1].translatedCodon == '3') {
        // Codon is not split
        left = intron.start - 2;
        right = intron.end + 2;
    } else {
        if (pairs[intron.start - 2].protein == '3'
                || pairs[intron.start - 2].translatedCodon == '3') {
            // Codon is split after the first nucleotide
            left = intron.start - 3;
            right = intron.end + 4;
            // Majority of the codon belongs to the right side
            double weight = kernel->getWeight(0);
            intron.rightScore += pairs[intron.end + 1].score(scoreMatrix) * weight;
            intron.rightWeightSum += weight;
        } else {
            // Codon is split after the second nucleotide
            left = intron.start - 4;
            right = intron.end + 3;
            // Majority of the codon belongs to the left side
            double weight = kernel->getWeight(0);
            intron.leftScore += pairs[intron.start - 1].score(scoreMatrix) * weight;
            intron.leftWeightSum += weight;
        }
    }

    scoreLeft(intron, left, windowWidth);
    scoreRight(intron, right, windowWidth);

    // Normalize alignments by their length, otherwise alignments
    // in short exons between introns are penalized
    if (multiply) {
        if (intron.leftWeightSum == 0 || intron.rightWeightSum == 0 ||
                intron.leftScore < 0 || intron.rightScore < 0) {
            intron.score = 0;
        } else {
            intron.score = (intron.leftScore / (intron.leftWeightSum)) *
                    (intron.rightScore / (intron.rightWeightSum));
            intron.score = sqrt(intron.score);
        }
    } else {
        intron.score = (intron.leftScore + intron.rightScore) /
                (intron.leftWeightSum + intron.rightWeightSum);
        if (intron.score < 0) {
            intron.score = 0;
        }
    }

    intron.scoreSet = true;
    return intron.score;
}

void Alignment::scoreLeft(Intron & intron, int start, int windowWidth) {
    for (int i = start; i > (start - windowWidth * 3); i -= 3) {
        // Check for end of local alignment
        if (i < 0 || pairs[i].type != 'e') {
            return;
        }
        double weight = kernel->getWeight((i - start) / 3);
        intron.leftWeightSum += weight;
        intron.leftScore += pairs[i].score(scoreMatrix) * weight;
    }
}

void Alignment::scoreRight(Intron & intron, int start, int windowWidth) {
    for (int i = start; i < (start + windowWidth * 3); i += 3) {
        // Check for end of local alignment
        if (i >= index || pairs[i].type != 'e') {
            return;
        }
        double weight = kernel->getWeight((i - start) / 3);
        intron.rightWeightSum += weight;
        intron.rightScore += pairs[i].score(scoreMatrix) * weight;
    }
}

Alignment::AlignedPair::AlignedPair(char tc, char n, char p) :
nucleotide(n),
translatedCodon(tc),
protein(p) {
    if (islower(n)) {
        this->type = 'i';
    } else {
        this->type = 'e';
    }
    // Assign an amino acid to a stop codon, so it is treated as a gap
    if (translatedCodon == '*') {
        translatedCodon = 'A';
    }

    // Fix strange J amino acid which is in fact S
    if (translatedCodon == 'J') {
        translatedCodon = 'S';
    }
}

double Alignment::AlignedPair::score(const ScoreMatrix * scoreMatrix) {
    return scoreMatrix->getScore(translatedCodon, protein);
}
