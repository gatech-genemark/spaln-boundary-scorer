#include <float.h>
#include <iostream>

#include "IntronStorage.h"

using namespace std;

IntronStorage::IntronStorage() {
    clear();
}

void IntronStorage::clear() {
    minScore = DBL_MAX;
    maxScore = 0;
    introns.clear();
}

void IntronStorage::storeIntron(string protein, string gene, int start, int end,
        char strand, string spliceSites, double score, int number, double leftExonScore,
        double rightExonScore) {

    Intron i(protein, gene, start, end, strand, spliceSites, score,
            number, leftExonScore, rightExonScore);
    introns.push_back(i);

    if (score > maxScore) {
        maxScore = score;
    }
    if (score < minScore) {
        minScore = score;
    }

}

void IntronStorage::normalizeScores() {
    normalizeScores(minScore, maxScore);
}

void IntronStorage::normalizeScores(double min, double max) {
    double range = max - min;
    if (range == 0) {
        return;
    }
    for (unsigned int i = 0; i < introns.size(); i++) {
        if (introns[i].score > max + 0.00001) {
            cerr << "Warning, supplied maximum score is lower than a score of "
                    "intron in " << introns[i].gene << "-" << introns[i].protein << endl;
            introns[i].score = max;
        }
        if (introns[i].score < min - 0.00001) {
            cerr << "Warning, supplied minimum score is higher than a score of "
                    "intron in " << introns[i].gene << "-" << introns[i].protein << endl;
            introns[i].score = min;
        }
        introns[i].score = (introns[i].score - min) / range;
    }
}

void IntronStorage::printIntrons(string output) {
    ofstream ofs(output.c_str(), std::ofstream::out | std::ofstream::app);
    for (unsigned int i = 0; i < introns.size(); i++) {
        ofs << introns[i].gene << "\tSpaln\tIntron\t";
        ofs << introns[i].start << "\t";
        ofs << introns[i].end << "\t";
        ofs << ".\t+\t.\tprot=" << introns[i].protein;
        ofs << "; intron_id=" << introns[i].number << ";";
        ofs << " splice_sites=" << introns[i].spliceSites << ";";
        ofs << " al_score=" << introns[i].score << ";";
        ofs << " LeScore=" << introns[i].leftExonScore << ";";
        ofs << " ReScore=" << introns[i].rightExonScore << ";" << endl;
    }
    ofs.close();
}
