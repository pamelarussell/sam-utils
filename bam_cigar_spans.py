
import argparse
import collections
import math
import re

import pybedtools

import matplotlib.pyplot as plt
import pysam
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
parser.add_argument('--span_res', action = 'store', dest = 'span_res', required = False, default = 100, help = 'Resolution for keeping track of cigar spans')
parser.add_argument('--log', action = 'store', dest = 'log', required = True, help = 'Log file')
args = parser.parse_args()
bam_file = args.bam
out_span_counts = args.out_counts
out_fig_prefix = args.out_fig_prefix
hist_label = args.hist_label
out_bam_prefix = args.out_bam_prefix
max_span = int(args.max_span)
chunk_size = int(args.chunk_size)
span_res = args.span_res 
log = args.log
logger = open(log, "w", buffering = 1)

# Determine the number of mapped reads in the bam file
logger.write("\nCounting mapped and unmapped reads in %s...\n" % bam_file)
n_mapped = int(pysam.view("-c", "-F", "4", bam_file))
n_unmapped = int(pysam.view("-c", "-f", "4", bam_file))
logger.write("There are %s mapped reads and %s unmapped reads.\n" % ("{:,}".format(n_mapped), "{:,}".format(n_unmapped)))

# Construct cigar length chunks for separate bam files
logger.write("\nDetermining cigar length chunks...\n")
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
logger.write("\nMaking bam writers...\n")
def chunk_to_bam(chunk):
    return "%s%s_%s.bam" % (out_bam_prefix, chunk[0], chunk[1])
bam_writers = {chunk: pysam.AlignmentFile(chunk_to_bam(chunk), "wb", header = header) for chunk in len_chunks}

# Keep track of cigar spans
# cigar_span_counts is a dict of dicts. Outer keys are reference names, values are dicts, inner dict maps span interval to count
# Keys of inner dict are the beginning of span interval, e.g., if the span intervals are size 10, a fragment with length 72 would
# count under the key 70
def span_to_key(span):
    return span - span % span_res
cigar_span_counts = {}
def add_cigar_span(ref, span):
    if ref not in cigar_span_counts:
        cigar_span_counts[ref] = collections.OrderedDict() 
        for s in range(max_span):
            cigar_span_counts[ref][s] = 0
    key = span_to_key(span)
    cigar_span_counts[ref][key] = 1 + cigar_span_counts[ref][key]
    
# Iterate through bam file and save cigar spans
# Write bam files of records with cigar span in each interval
logger.write("\nIterating through bam file and getting cigar spans...\n")
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
    if i % 1000000 == 0:
        logger.write("Finished %s records\n" % "{:,}".format(i))
logger.write("Finished iterating through bam file.\n")

# Close the bam files
logger.write("\nClosing bam reader...\n")
bam_reader.close()
logger.write("\nClosing the bam writers...\n")
for bam_writer in bam_writers.values():
    bam_writer.close()

# Index the chunk bam files
logger.write("\nIndexing the chunk bam files\n")
for chunk in len_chunks:
    logger.write("Chunk $s/%s\n" % (chunk, len_chunks))
    pysam.index(chunk_to_bam(chunk))
    
# Write bedgraph coverage files for each chunk bam file for easier display in IGV
logger.write("\nWriting bedgraph coverage files\n")
def chunk_to_bedgraph(chunk):
    return "%s%s_%s.bedgraph" % (out_bam_prefix, chunk[0], chunk[1])
for chunk in len_chunks:
    logger.write("Chunk $s/%s\n" % (chunk, len_chunks))
    bedtool = pybedtools.BedTool(chunk_to_bam(chunk))
    cov = bedtool.genome_coverage(bg = True, split = True)
    cov.saveas(chunk_to_bedgraph(chunk))
    
# Write the span counts to a table
if out_span_counts is not None:
    logger.write("\nWriting span counts to file:\n%s\n" % out_span_counts)
    with open(out_span_counts, 'w') as w:
        w.write("ref\tspan\tnum_records\n")
        for ref, d in cigar_span_counts.items():
            for span, count in d.items():
                if count > 0:
                    w.write("%s\t%s\t%s\n" % (ref, span, count))
        
# Save histograms of cigar spans for each reference sequence
def to_log(x):
    if x == 0:
        return 0
    if x > 0:
        return math.log10(x)
    else:
        raise ValueError("Invalid arg to log: %s" %x)
if out_fig_prefix is not None:
    logger.write("\n")
    for ref in cigar_span_counts.keys():
        out_fig = "%s%s.pdf" % (out_fig_prefix, re.sub("[| .]", r'_', ref))       
        logger.write("Writing histogram of cigar spans to file:\n%s\n" % out_fig)
        plt_data = [to_log(x) for x in cigar_span_counts[ref].values()]
        plt.bar(range(max_span), plt_data)
        plt_title = ref
        if hist_label is not None:
            plt_title = "%s -> %s" % (hist_label, plt_title)
        plt.title(plt_title)
        plt.xlabel("Cigar span")
        plt.ylabel("Number of reads (log10)")
        plt.savefig(out_fig)
    
    
logger.write("\nAll done.\n\n")
logger.close()    

