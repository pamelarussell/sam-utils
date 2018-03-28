# sam-utils: Utilities for analyzing SAM/BAM files

## bam_cigar_spans.py

This script reads a bam file and determines the span along the reference sequence covered by each record. A table is written that gives the number of records with each span for each reference sequence. For each reference sequence, a histogram of read spans is also saved. The output files are optional and are not generated if output file paths are not specified.

To run:

`python bam_cigar_spans.py --help`
