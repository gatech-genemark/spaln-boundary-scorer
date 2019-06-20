#ifndef ALIGNMENT_H
#define ALIGNMENT_H

#include <fstream>
#include <vector>
#include <string>
#include "ScoreMatrix.h"
#include "IntronStorage.h"
#include "Kernel.h"

using namespace std;

/// Class for parsing a single gene-protein alignment

class Alignment {
public:
    Alignment();
    ~Alignment();
    /**
     * Parse a single gene-protein alignment.
     * The function checks if the general structure of the alignment
     * is OK but it does not check the validity of every single base/protein.
     *
     * @param fstream File stream starting at the position of
     *                the alignment start.
     */
    int parse(ifstream & inputStream, string headerLine);
    /**
     * @return Name of the aligned gene
     */
    string getGene();
    /**
     * @return Name of the aligned protein
     */
    string getProtein();
    /**
     * @return Total alignment length in nucleotides
     */
    int getLength();
    /**
     * Print the whole alignment
     */
    void print(ostream & os);
    /**
     * Store introns found in this alignment to the IntronStorage object
     * which contains all introns
     * @param storeage The intron storage
     */
    void storeIntrons(IntronStorage & storage);
    /**
     * Score all introns in the alignment
     * @param windowWidth Number of amino acids scored in the upstream and
     *                    downstream regions
     * @param multiply    Whether to multiply or sum scores from upstream and
     *                    downstream regions
     * @param scoreMatrix Scoring matrix used for scoring amino acids. If no
     *                    matrix is given, score is determined based on the
     *                    quality indicator in the input alignment file
     */
    void scoreIntrons(int windowWidth, bool multiply,
            const ScoreMatrix * scoreMatrix, Kernel * kernel);
    /**
     * @return True if alignment contains any introns
     */
    bool hasIntrons();
private:
    /// Starting position of the alignment in dna
    int dnaStart;
    /// Starting position of the alignment in protein
    int proteinStart;
    int i;
    string gene;
    string protein;
    /**
     * Clear the object for new alignment pair
     */
    void clear();

    /// Overall alignment lentgh
    unsigned int blockLength;
    /// Overall position in alignment
    int index;
    int realPositionCounter;

    /// Single nucleotide-protein pair
    struct AlignedPair {
        /**
         * Save pair and determine exon/intron
         */
        AlignedPair(char tc, char n, char p, bool insideIntron);
        /**
         * Return amino acid score
         * @return AA score
         */
        double score(const ScoreMatrix * scoreMatrix);
        char nucleotide;
        /**
         * The proteins are saved as follows: 1P3 or ppp
         * Where P is the protein, numbers 1 and 3 fill
         * the space created by 3 nucleotide to 1 protein mapping.
         * Lowercase letters are used at splice sites.
         */
        char translatedCodon;
        char protein;
        /**
         * Type of DNA base (based on alignment)
         * 'i' for intron
         * 'e' for exon
         * ' ' for unknown
         */
        char type;
        /// Position of a nucleotide in the alignment relative to a gene start
        int realPosition;
    };

    struct Exon {
        Exon(int start);
        unsigned int start, end;
        double score;
        bool scoreSet;
    };

    /**
     * Intron structure for the purposes of the alignment only. Final set of
     * introns from all alignments with correct start and end positions
     * relative to gene start is stored in the IntronStorage class
     */
    struct Intron {
        Intron();
        unsigned int start, end;
        double score;
        bool scoreSet;
        char donor[2];
        char acceptor[2];
        bool complete;
        double leftScore, rightScore;
        /// Flag indicating that a gap, or aligned protein, was detected inside intron
        bool gap;
        Exon * leftExon;
        Exon * rightExon;
    };

    /// Initial size of alignment vector
    static const int N = 3000;
    /// Array containing the actual alignment
    vector<AlignedPair> pairs;


    static const int BLOCK_LENGTH = 100;
    static const int BLOCK_ITEMS_CNT = 3;
    static const int BLOCK_OFFSET = 9;
    /**
     *  Parse gene and protein name from header line
     */
    int parseHeader(string headerLine);
    /**
     *  Parse individual block of lines containing the alignment and its properties
     */
    void parseBlock(const vector<string>& lines);

    /**
     * Check if the given character is an amino acid or a gap
     */
    bool gapOrAA(char a);
    /**
     * Assign phases to all positions
     */
    void assignCodonPhases();
    void printLineError();

    bool insideIntron;
    // Flag indicating that donor position of an intron is being read
    bool donorFlag;
    /**
     * Detect and save introns.
     * The function also retrieves information associated with the intron
     * such as its start/end position and donor/acceptor site.
     * @param pair
     */
    void checkForIntron(AlignedPair & pair);

    void checkForStart(AlignedPair & pair);
    vector<Intron> introns;
    vector<Exon*> exons;
    const ScoreMatrix * scoreMatrix;
    Kernel * kernel;

    /**
     * Determine score of a single intron using exon alignment in the
     * upstream and downstream region
     */
    double scoreIntron(Intron & intron, int windowWidth, bool multiply);
    /// Alignment score of amino acids in exon before intron start
    void scoreLeft(Intron & intron, int start, int windowWidth);
    /// Alignment score of amino acids in exon after intron end
    void scoreRight(Intron & intron, int start, int windowWidth);
    void scoreExon(Exon * exon);
};


#endif /* ALIGNMENT_H */


