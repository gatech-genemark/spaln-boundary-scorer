#ifndef SCORE_MATRIX_H
#define SCORE_MATRIX_H

#include <string>
#include <fstream>
#include <map>
#include <vector>

#define UNKNOWN_SCORE -4

using namespace std;

// Class for loading and accessing substitution scoring matrices

class ScoreMatrix {
public:
    /// Load scoring matrix from a file in csv format
    bool loadFromFile(string filename);
    double getScore(char a, char b) const;
    double getMaxScore() const;
    void print() const;
private:
    map<char, map<char, double> > matrix;
    ifstream inputStream;
    int size;
    vector<char> columnHeaders;
    bool readColumnHeaders();
    bool readRow();
    // Replace common csv separators by whitespaces to make further
    // parsing easier
    void processLine(string & line);
    void computeMaxScore();
    double maxScore;
};

#endif /* SCORE_MATRIX_H */
