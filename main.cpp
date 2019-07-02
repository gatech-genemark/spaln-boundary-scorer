#include "Parser.h"
#include "ScoreMatrix.h"
#include "Kernel.h"

#include <iostream>
#include <string>
#include <unistd.h>
#include <stdlib.h>

using namespace std;

#define DEFAULT_WINDOW_WIDTH 10
#define DEFAULT_KERNEL "triangular"
#define DEFAULT_EXON_SCORE 25

void printUsage(char * name) {
    cout << "Usage: " << name << " < input -o output_file "
            "[-w integer] [-am] [-s matrix_file] [-k kernel] [-e min_exon_score]" << endl;
    cout << "Options:" << endl;
    cout << "   -w Width of a scoring window around introns. Default = " <<
            DEFAULT_WINDOW_WIDTH << endl;
    cout << "   -s Use specified scoring matrix for intron scoring. Wi-\n"
            "      thout this option, amino acid scores are determined \n"
            "      based on the quality indicator in the input alignment \n"
            "      file." << endl;
    cout << "   -k Specify type of weighting kernel used. Available opti-\n"
            "      ons are \"triangular\", \"box\", \"parabolic\" and \n"
            "      \"triweight\". Triangular kernel is the default option." << endl;
    cout << "   -e Minimum exon score. Exons with lower scores (and int-\n"
            "      rons bordering such exons; start and stops inside the \n"
            "      exons) are not printed. Default = " <<
            DEFAULT_EXON_SCORE << endl;
}

int main(int argc, char** argv) {
    int opt;
    int windowWidth = DEFAULT_WINDOW_WIDTH;
    string output;
    string matrixFile = "";
    string kernelType = DEFAULT_KERNEL;
    double minExonScore = DEFAULT_EXON_SCORE;

    while ((opt = getopt(argc, argv, "o:w:s:k:e:")) != EOF) {
        switch (opt) {
            case 'o':
                output = optarg;
                break;
            case 'w':
                windowWidth = atoi(optarg);
                break;
            case 's':
                matrixFile = optarg;
                break;
            case 'k':
                kernelType = optarg;
                break;
            case 'e':
                minExonScore = atof(optarg);
                break;
            case '?':
                printUsage(argv[0]);
                return 1;
            default:
                return 1;
        }

    }

    if (output.size() == 0) {
        cerr << "error: Output file not specified" << endl;
        printUsage(argv[0]);
        return 1;
    }

    if (kernelType != "triangular" && kernelType != "box" &&
            kernelType != "parabolic" && kernelType != "triweight") {
        cerr << "error: Invalid kernel. Valid options are \"box\","
                "\"triangular\", \"parabolic\" and \"triweight\" kernels." << endl;
        printUsage(argv[0]);
        return 1;
    }

    ScoreMatrix * scoreMatrix = NULL;
    if (!matrixFile.empty()) {
        scoreMatrix = new ScoreMatrix();
        if (!scoreMatrix->loadFromFile(matrixFile)) {
            cerr << "error: Could not load scoring matrix" << endl;
            printUsage(argv[0]);
            return 1;
        }
    }

    Kernel * kernel;
    if (kernelType == "triangular") {
        kernel = new TriangularKernel();
    } else if (kernelType == "box") {
        kernel = new BoxKernel();
    } else if (kernelType == "parabolic") {
        kernel = new ParabolicKernel();
    } else if (kernelType == "triweight") {
        kernel = new TriweightKernel();
    }

    Parser fileParser;
    fileParser.setWindowLegth(windowWidth);
    fileParser.setScoringMatrix(scoreMatrix);
    fileParser.setKernel(kernel);
    fileParser.setMinExonScore(minExonScore);

    int result = fileParser.parse(output);

    delete scoreMatrix;
    delete kernel;
    return result;
}

