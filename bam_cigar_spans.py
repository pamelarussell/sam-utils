"""
bam_cigar_spans.py

This script reads a bam file and determines the span along the reference sequence
covered by each record. A table is written that gives the number of records with
each span for each reference sequence. For each reference sequence, a histogram
of read spans is also saved. The output files are optional and are not generated
if output file paths are not specified.
"""

import argparse
import re

import pysam

import matplotlib.pyplot as plt
import pandas as pd
from sam_utils import cigar_span

# Parse the command line
parser = argparse.ArgumentParser()
parser.add_argument('--bam', action = 'store', dest = 'bam', required = True, help = 'Bam file')
parser.add_argument('--out_counts', action = 'store', dest = 'out_counts', required = False, help = 'Output table of cigar span counts')
parser.add_argument('--out_hist_prefix', action = 'store', dest = 'out_fig_prefix', required = False, help = 'Prefix for output histograms (one per reference sequence)')
parser.add_argument('--extra_hist_label', action = 'store', dest = 'hist_label', required = False, help = 'Additional label to prepend to histogram titles')
args = parser.parse_args()
bam_file = args.bam
out_span_counts = args.out_counts
out_fig_prefix = args.out_fig_prefix
hist_label = args.hist_label

# Open the bam file and get an iterator over records
print("Opening bam file:\n%s" % bam_file)
bam_reader = pysam.AlignmentFile(bam_file, "rb")
bam_iter = bam_reader.fetch()

# Determine the number of mapped reads in the bam file
n_mapped = int(pysam.view("-c", "-F", "4", bam_file))
n_unmapped = int(pysam.view("-c", "-f", "4", bam_file))
print("\nThere are %s mapped reads and %s unmapped reads." % (n_mapped, n_unmapped))

# Keep track of cigar spans
cigar_spans = pd.DataFrame(columns = ["ref", "span"], index = range(n_mapped))

# For each read, record the cigar span in the map corresponding to its reference sequence
print("\nIterating through bam file and getting cigar spans...")
i = 0
for rec in bam_iter:
    ref = rec.reference_name
    span = cigar_span(rec.cigartuples)
    cigar_spans.iloc[i] = [ref, span]
    i = i + 1
    if i % 10000 == 0:
        print("Finished %s records" % i)

# Close the bam file
bam_reader.close()

# Write the span counts to a table
if out_span_counts is not None:
    print("\nWriting span counts to file:\n%s" % out_span_counts)
    cigar_spans.\
        groupby(["ref", "span"]).\
        aggregate(len).\
        reset_index().\
        rename(columns = {0: "num_records"}).\
        to_csv(out_span_counts, sep = "\t", index = False)
        
# Save histograms of cigar spans for each reference sequence
if out_fig_prefix is not None:
    print("")
    for ref in cigar_spans.dropna().ref.unique():
        out_fig = "%s%s.pdf" % (out_fig_prefix, re.sub("[| .]", r'_', ref))
        fig_data = cigar_spans.loc[cigar_spans["ref"] == ref].as_matrix(columns = ["span"])
        print("Writing histogram of cigar spans for %s reads to file:\n%s" % (fig_data.shape[0], out_fig))
        n, bins, patches = plt.hist(fig_data, bins = 500, histtype = "stepfilled", cumulative = False)
        plt_title = ref
        if hist_label is not None:
            plt_title = "%s -> %s" % (hist_label, plt_title)
        plt.title(plt_title)
        plt.xlabel("Cigar span")
        plt.ylabel("Number of reads")
        plt.savefig(out_fig)
    
    
print("\nAll done.\n")
    

