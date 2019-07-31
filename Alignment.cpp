#include "Alignment.h"
#include "Parser.h"
#include <string>
#include <iostream>
#include <fstream>
#include <sstream>
#include <utility>
#include <cmath>
#include <ctype.h>
#include <algorithm>

using namespace std;

Alignment::Alignment() {
    pairs.reserve(N);
    start = NULL;
    stop = NULL;
}

Alignment::~Alignment() {
    clear();
}

void Alignment::clear() {
    index = 0;
    insideIntron = false;
    donorFlag = false;
    introns.clear();
    for (unsigned int i = 0; i < exons.size(); i++) {
        delete exons[i];
    }
    exons.clear();
    dnaStart = 0;
    if (start != NULL) {
        delete start;
    }
    start = NULL;
    if (stop != NULL) {
        delete stop;
    }
    stop = NULL;
}

int Alignment::parse(istream& inputStream, string headerLine, bool forward) {
    clear();
    this->forward = forward;
    vector<string> blockLines(BLOCK_ITEMS_CNT);
    string line;

    // Read header
    int status = parseHeader(headerLine);
    if (status != READ_SUCCESS) {
        cerr << "error: Invalid alignment header " << endl;
        return status;
    }
  //  cout << gene << " " << protein << endl;
    // Get to the alignment
    getline(inputStream, line);
    while (line.substr(0, 9) != "ALIGNMENT") {
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
                dnaStart = atoi(line.substr(0, BLOCK_OFFSET).c_str());
                if (dnaStart <= 0) {
                    cerr << "error: Could not read dna alignment start position" << endl;
                    return FORMAT_FAIL;
                }
                realPositionCounter = dnaStart;
            } else if (i == 2) {
                proteinStart = atoi(line.substr(0, BLOCK_OFFSET).c_str());
                if (proteinStart <= 0) {
                    cerr << "error: Could not read dna alignment start position" << endl;
                    return FORMAT_FAIL;
                }
            }
        } else {
            printLineError();
            return FORMAT_FAIL;
        }
    }
    // Any symbol in DNA other than space before pipe signifies the true end of
    // alignment Looking for the first whitespace in DNA does not work as there
    // might be gaps inside introns
    int pipePosition = blockLines[1].find("|");
    blockLength = blockLines[1].find_last_not_of(" ", pipePosition - 1) - BLOCK_OFFSET  + 1;
    for (i = 0; i < BLOCK_ITEMS_CNT; i++) {
        blockLines[i] = blockLines[i].substr(BLOCK_OFFSET, blockLength);
        if (blockLines[i].size() != blockLength) {
            printLineError();
            return FORMAT_FAIL;
        }
    }

    replace(blockLines[1].begin(), blockLines[1].end(), ' ', '-');
    parseBlock(blockLines);


    return READ_SUCCESS;
}

int Alignment::parseHeader(string headerLine) {
    bool geneParsed = false;
    for (unsigned int i = 0; i < headerLine.size(); i++) {
        if (headerLine[i] == '>' || headerLine[i] == '<') {
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
        AlignedPair pair(lines[0][i], lines[1][i], lines[2][i], insideIntron);

        checkForIntron(pair);
        checkForStart(pair);
        checkForStop(pair);

        if (pair.nucleotide != '-' && pair.nucleotide != ' ') {
            if (forward) {
                pair.realPosition = realPositionCounter++;
            } else {
                pair.realPosition = realPositionCounter--;
            }
        } else {
            // Don't increment the counter
            if (forward) {
                pair.realPosition = realPositionCounter - 1;
            } else {
                pair.realPosition = realPositionCounter + 1;
            }
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
            } else if (i == blockLength - 1) {
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
            } else if (i == blockLength - 1) {
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
    // Some alignments start with a gap or intron, do not create initial exon
    // in such cases
    int alignmentPosition;
    if (forward) {
        alignmentPosition = realPositionCounter - dnaStart + 1;
    } else {
        alignmentPosition = dnaStart - realPositionCounter + 1;
    }
    if (alignmentPosition == 1 && pair.type == 'e' && pair.nucleotide != '-') {
        exons.push_back(new Exon(index));
    }

    if (index == (int) blockLength - 1 && pair.type == 'e') {
        exons.back()->end = index;
    }

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
        if (index != 0) {
            exons.back()->end = index - 1;
            i.leftExon = exons.back();
        }
        introns.push_back(i);
        insideIntron = true;
        donorFlag = true;
    } else if (insideIntron && pair.type != 'i') { // intron end
        introns.back().end = index - 1;
        introns.back().acceptor[0] = pairs[index - 2].nucleotide;
        introns.back().acceptor[1] = pairs[index - 1].nucleotide;
        insideIntron = false;
        if (introns.back().start != 0 && introns.back().gap == false) {
            introns.back().complete = true;
        }
        exons.push_back(new Exon(index));
        introns.back().rightExon = exons.back();
    } else if (insideIntron && pair.nucleotide == '-') {
        // Gap (AA aligned) inside introns, do not report these introns
        introns.back().gap = true;
    }
}

void Alignment::checkForStart(AlignedPair& pair) {
    //int alignmentPosition = realPositionCounter - dnaStart + 1;
    if (index == 2) {
        string codon = "";
        codon += pair.nucleotide;
        for (int i = 1; i < 3; i++) {
            codon = pairs[index - i].nucleotide + codon;
        }
        if (codon == "ATG") {
            // Check if protein alignment starts with its first M
            if (proteinStart == 1 && pairs[index - 1].protein == 'M') {
                start = new Codon(index - 2, exons.back());
            }
        }
    }
}

void Alignment::checkForStop(AlignedPair& pair) {
    if (index == (int) blockLength - 1 && pairs[index - 3].type == 'e' &&
            (pair.nucleotide == 'a' || pair.nucleotide == 'g')) {
        string codon = "";
        codon += pair.nucleotide;
        for (int i = 1; i < 3; i++) {
            codon = pairs[index - i].nucleotide + codon;
        }

        if (codon == "taa" || codon == "tag" || codon == "tga") {
            stop = new Codon(index - 2, exons.back());
        }
    }
}

void Alignment::printHints(string output, double minExonScore) {
    ofstream ofs(output.c_str(), std::ofstream::out | std::ofstream::app);
    char strand;
    if (forward) {
        strand = '+';
    } else {
        strand = '-';
    }

    for (unsigned int i = 0; i < introns.size(); i++) {
        if (!introns[i].complete || introns[i].leftExon->score < minExonScore ||
                introns[i].rightExon->score < minExonScore) {
            continue;
        }

        string spliceSites(introns[i].donor, 2);
        spliceSites.append("_");
        spliceSites.append(introns[i].acceptor, 2);

        ofs << gene << "\tSpaln_Parser\tIntron\t";
        if (forward) {
            ofs << pairs[introns[i].start].realPosition << "\t";
            ofs << pairs[introns[i].end].realPosition << "\t";
        } else {
            ofs << pairs[introns[i].end].realPosition << "\t";
            ofs << pairs[introns[i].start].realPosition << "\t";
        }
        ofs << ".\t" << strand << "\t.\tprot=" << protein;
        ofs << "; intron_id=" << i + 1 << ";";
        ofs << " splice_sites=" << spliceSites << ";";
        ofs << " al_score=" << introns[i].score << ";";
        ofs << " LeScore=" << introns[i].leftExon->score << ";";
        ofs << " ReScore=" << introns[i].rightExon->score << ";\n";
    }

    if (start != NULL && start->exon->score >= minExonScore) {
        ofs << gene << "\tSpaln_Parser\tstart_codon\t";
        if (forward) {
            ofs << pairs[start->position].realPosition << "\t";
            ofs << pairs[start->position + 2].realPosition  << "\t";
        } else {
            ofs << pairs[start->position + 2].realPosition  << "\t";
            ofs << pairs[start->position].realPosition << "\t";
        }
        ofs << ".\t" << strand << "\t.\tprot=" << protein << ";";
        ofs << " score=" << start->score << ";";
        ofs << " eScore=" << start->exon->score << ";\n";
    }

    for (unsigned int i = 0; i < exons.size(); i++) {
        if (exons[i]->score < minExonScore) {
            continue;
        }
        ofs << gene << "\tSpaln_Parser\tCDS\t";
        if (forward) {
            ofs << pairs[exons[i]->start].realPosition << "\t";
            ofs << pairs[exons[i]->end].realPosition << "\t";
        } else {
            ofs << pairs[exons[i]->end].realPosition << "\t";
            ofs << pairs[exons[i]->start].realPosition << "\t";
        }
        ofs << ".\t" << strand << "\t.\tprot=" << protein;
        ofs << "; exon_id=" << i + 1 << ";";
        ofs << " score=" << exons[i]->score << ";\n";
    }

    if (stop != NULL && stop->exon->score >= minExonScore) {
        ofs << gene << "\tSpaln_Parser\tstop_codon\t";
        if (forward) {
            ofs << pairs[stop->position].realPosition << "\t";
            ofs << pairs[stop->position + 2].realPosition  << "\t";
        } else {
            ofs << pairs[stop->position + 2].realPosition  << "\t";
            ofs << pairs[stop->position].realPosition << "\t";
        }
        ofs << ".\t" << strand << "\t.\tprot=" << protein << ";";
        ofs << " score=" << stop->score << ";";
        ofs << " eScore=" << stop->exon->score << ";\n";
    }

    ofs.close();
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
    gap = false;
}

Alignment::Exon::Exon(int start) {
    this->start = start;
    scoreSet = false;
}

Alignment::Codon::Codon(int position, Exon * exon) {
    this->position = position;
    this->exon = exon;
}

void Alignment::scoreHints(int windowWidth,
        const ScoreMatrix * scoreMatrix, Kernel * kernel) {
    this->scoreMatrix = scoreMatrix;
    this->kernel = kernel;
    this->kernel->setWidth(windowWidth);

    for (unsigned int i = 0; i < exons.size(); i++) {
        scoreExon(exons[i]);
    }

    if (start != NULL) {
        start->score = 0;
        scoreStart(start, windowWidth);
        start->score /=  kernel->weightSum();
        start->score /= scoreMatrix->getMaxScore();
    }

    if (stop != NULL) {
        stop->score = 0;
        scoreStop(stop, windowWidth);
        stop->score /=  kernel->weightSum();
        stop->score /= scoreMatrix->getMaxScore();
    }

    for (unsigned int i = 0; i < introns.size(); i++) {
        if (introns[i].complete && !introns[i].scoreSet) {
            scoreIntron(introns[i], windowWidth);
        }
    }

}

double Alignment::scoreIntron(Intron& intron, int windowWidth) {
    intron.leftScore = intron.rightScore = 0;
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
            right = intron.end + 1;
        } else {
            // Codon is split after the second nucleotide
            left = intron.start - 1;
            right = intron.end + 3;
        }
    }

    scoreLeft(intron, left, windowWidth);
    scoreRight(intron, right, windowWidth);
    double weightSum = kernel->weightSum();

    // Normalize alignments by the area under kernel

    if (intron.leftScore <= 0 || intron.rightScore <= 0) {
        intron.score = 0;
    } else {
        intron.score = (intron.leftScore / weightSum) *
                (intron.rightScore / weightSum);
        intron.score = sqrt(intron.score);
    }

    intron.score /= scoreMatrix->getMaxScore();
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
        intron.rightScore += pairs[i].score(scoreMatrix) * weight;
    }
}

void Alignment::scoreStart(Codon* start, int windowWidth) {
    for (int i = start->position + 1; i < (start->position + 1 + windowWidth * 3); i += 3) {
        // Check for end of local alignment
        if (i >= index || pairs[i].type != 'e') {
            return;
        }
        double weight = kernel->getWeight((i - start->position - 1) / 3);
        start->score += pairs[i].score(scoreMatrix) * weight;
    }
}

void Alignment::scoreStop(Codon* stop, int windowWidth) {
    for (int i = stop->position - 2; i > (stop->position - 2 - windowWidth * 3); i -= 3) {
        // Check for end of local alignment
        if (i < 0 || pairs[i].type != 'e') {
            return;
        }
        double weight = kernel->getWeight((i - stop->position + 2) / 3);
        stop->score += pairs[i].score(scoreMatrix) * weight;
    }
}


void Alignment::scoreExon(Exon* exon) {
    int i = exon->start;
    exon->score = 0;
    while (i <= exon->end) {
        if (gapOrAA(pairs[i].protein)) {
            exon->score += pairs[i].score(scoreMatrix);
        }
        i++;
    }
}

Alignment::AlignedPair::AlignedPair(char tc, char n, char p, bool insideIntron) :
nucleotide(n),
translatedCodon(tc),
protein(p) {
    if (islower(n) || (insideIntron && n == '-')) {
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
