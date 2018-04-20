
import argparse
import re

import pysam

import matplotlib.pyplot as plt
import pandas as pd
import pybedtools
from sam_utils import cigar_span


# Parse the command line
parser = argparse.ArgumentParser(description = """
bam_cigar_spans.py
 
This script reads a bam file and determines the span along the reference sequence
covered by each record. It creates three sets of outputs:

1. A table giving the number of records with each span for each reference sequence
2. A histogram of read spans for each reference sequence
3. Bam files consisting only of reads whose span is in a particular length interval
""")
parser.add_argument('--bam', action = 'store', dest = 'bam', required = True, help = 'Bam file')
parser.add_argument('--out_counts', action = 'store', dest = 'out_counts', required = False, help = 'Output table of cigar span counts')
parser.add_argument('--out_fig_prefix', action = 'store', dest = 'out_fig_prefix', required = False, help = 'Prefix for output histograms (one per reference sequence)')
parser.add_argument('--extra_hist_label', action = 'store', dest = 'hist_label', required = False, help = 'Additional label to prepend to histogram titles')
parser.add_argument('--out_bam_prefix', action = 'store', dest = 'out_bam_prefix', required = True, help = 'Prefix for output bam files (one file per read span chunk)')
parser.add_argument('--max_read_span_out', action = 'store', dest = 'max_span', required = True, help = 'Max read span to count / write to chunk bam files')
parser.add_argument('--chunk_size_out', action = 'store', dest = 'chunk_size', required = True, help = 'Size of chunks (length interval)')
args = parser.parse_args()
bam_file = args.bam
out_span_counts = args.out_counts
out_fig_prefix = args.out_fig_prefix
hist_label = args.hist_label
out_bam_prefix = args.out_bam_prefix
max_span = int(args.max_span)
chunk_size = int(args.chunk_size)

# Determine the number of mapped reads in the bam file
n_mapped = int(pysam.view("-c", "-F", "4", bam_file))
n_unmapped = int(pysam.view("-c", "-f", "4", bam_file))
print("\nThere are %s mapped reads and %s unmapped reads." % (n_mapped, n_unmapped))

# Construct cigar length chunks for separate bam files
len_chunks = []
begin = 0
end = chunk_size
while end <= max_span:
    len_chunks.append((begin, end))
    begin = end
    end = begin + chunk_size
if end < max_span:
    len_chunks.append((end, max_span))
# Map of cigar length to length chunk
len_to_chunk = {l: t for t in len_chunks for l in range(t[0], t[1])}

# Open the bam file and get the header
bam_reader = pysam.AlignmentFile(bam_file, "rb")
header = bam_reader.header

# Make a bam file writer for each length range/chunk
def chunk_to_bam(chunk):
    return "%s%s_%s.bam" % (out_bam_prefix, chunk[0], chunk[1])
bam_writers = {chunk: pysam.AlignmentFile(chunk_to_bam(chunk), "wb", header = header) for chunk in len_chunks}

# Keep track of cigar spans
cigar_span_counts = {}
def add_cigar_span(ref, span):
    if ref not in cigar_span_counts:
        cigar_span_counts[ref] = {s: 0 for s in range(max_span)}
    cigar_span_counts[ref][span] = 1 + cigar_span_counts[ref][span]
    
# Iterate through bam file and save cigar spans
# Write bam files of records with cigar span in each interval
print("\nIterating through bam file and getting cigar spans...")
bam_iter = bam_reader.fetch()
i = 0
for rec in bam_iter:
    ref = rec.reference_name
    span = cigar_span(rec.cigartuples)
    if span > max_span:
        continue
    # Append the record to the appropriate bam file
    bam_writers[len_to_chunk[span]].write(rec)
    add_cigar_span(ref, span)
    i = i + 1
    if i % 10000 == 0:
        print("Finished %s records" % i)

# Close the bam files
bam_reader.close()
for bam_writer in bam_writers.values():
    bam_writer.close()

# Index the chunk bam files
for chunk in len_chunks:
    pysam.index(chunk_to_bam(chunk))
    
# Write bedgraph coverage files for each chunk bam file for easier display in IGV
def chunk_to_bedgraph(chunk):
    return "%s%s_%s.bedgraph" % (out_bam_prefix, chunk[0], chunk[1])
for chunk in len_chunks:
    bedtool = pybedtools.BedTool(chunk_to_bam(chunk))
    cov = bedtool.genome_coverage(bg = True, split = True)
    cov.saveas(chunk_to_bedgraph(chunk))
    
# Write the span counts to a table
if out_span_counts is not None:
    print("\nWriting span counts to file:\n%s" % out_span_counts)
    with open(out_span_counts, 'w') as w:
        w.write("ref\tspan\tnum_records\n")
        for ref, d in cigar_span_counts.items():
            for span, count in d.items():
                if count > 0:
                    w.write("%s\t%s\t%s\n" % (ref, span, count))
        
# Save histograms of cigar spans for each reference sequence
if out_fig_prefix is not None:
    raise RuntimeError("Histogram not implemented")
# cigar_spans was a data frame with columns ref, span and a row for each mapped read
#     print("")
#     for ref in cigar_spans.dropna().ref.unique():
#         out_fig = "%s%s.pdf" % (out_fig_prefix, re.sub("[| .]", r'_', ref))
#         fig_data = cigar_spans.loc[cigar_spans["ref"] == ref].as_matrix(columns = ["span"])
#         print("Writing histogram of cigar spans for %s reads to file:\n%s" % (fig_data.shape[0], out_fig))
#         n, bins, patches = plt.hist(fig_data, bins = 500, histtype = "stepfilled", cumulative = False, log = True)
#         plt_title = ref
#         if hist_label is not None:
#             plt_title = "%s -> %s" % (hist_label, plt_title)
#         plt.title(plt_title)
#         plt.xlabel("Cigar span")
#         plt.ylabel("Number of reads")
#         plt.savefig(out_fig)
    
    
print("\nAll done.\n")
    

