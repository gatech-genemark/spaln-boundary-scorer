# Description

Spaln boundary scorer parses introns, starts, stops and exons from Spaln's
alignment output and scores them. Introns, starts and stops are scored based
on local alignment quality around their boundaries. Detailed description of
how the scores are computed is available in:

## Installation and usage

To install, clone this repository and run `make` in the root folder.

To run, use the following command:

    spaln_boundary_scorer < spaln_input -o output_file -s matrix_file [-w integer] [-k kernel] [-e min_exon_score] [-r]

Input details:

* The program can parse multiple separate alignments saved in the same input.
* The input is read from stdin.
* Each input alignment is assumed to be on a single line (number of characters
per line, controlled by `-l` option in Spaln, is larger than the alignment
length).

Available options are:

```
   -o Where to save output file
   -s Path to amino acid scoring matrix
   -w Width of a scoring window around introns. Default = 10
   -k Specify type of weighting kernel used. Available opti-
      ons are "triangular", "box", "parabolic" and 
      "triweight". Triangular kernel is the default option.
   -e Minimum exon score. Exons with lower scores (and int-
      rons bordering such exons; start and stops inside the 
      exons) are not printed. Default = 25
   -r Process alignments on the reverse DNA strand (which are
      ignored by default)
```

## Tests

Unit tests are located in the `test` folder. To compile a test binary, run
`make test` in the root folder. Subsequently, run the test binary to evaluate
the tests:

    make test
    test/t_spaln_boundary_scorer
