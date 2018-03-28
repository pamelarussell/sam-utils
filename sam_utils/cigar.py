
class _CigarOperation():
    def __init__(self, op, bam, consumes_query, consumes_ref):
        """ Instantiate a cigar element.
        
        Args:
            op: Operation (character e.g. 'M')
            bam: Numeric identifier e.g. 0 for 'M'
            consumes_query: Whether the operation consumes the query sequence (boolean)
            consumes_ref: Whether the operation consumes the reference sequence (boolean)
        """
        self.op = op
        self.bam = bam
        self.consumes_query = consumes_query
        self.consumes_ref = consumes_ref

# The cigar operations
cigar_align_match = _CigarOperation('M', 0, True, True)
cigar_insertion = _CigarOperation('I', 1, True, False)
cigar_deletion = _CigarOperation('D', 2, False, True)
cigar_skip = _CigarOperation('N', 3, False, True)
cigar_soft_clip = _CigarOperation('S', 4, True, False)
cigar_hard_clip = _CigarOperation('H', 5, False, False)
cigar_pad = _CigarOperation('P', 6, False, False)
cigar_seq_match = _CigarOperation('=', 7, True, True)
cigar_seq_mismatch = _CigarOperation('X', 8, True, True)

# List of cigar operations
cigar_ops = [cigar_align_match,
                  cigar_insertion,
                  cigar_deletion,
                  cigar_skip,
                  cigar_soft_clip,
                  cigar_hard_clip,
                  cigar_pad,
                  cigar_seq_match,
                  cigar_seq_mismatch]

# Map of bam number to cigar operation
bam_to_cigar_element = {elt.bam: elt for elt in cigar_ops}

# Map of character to cigar operation
op_to_cigar_element = {elt.op: elt for elt in cigar_ops}

# Map of bam number to whether the operation consumes reference sequence
bam_to_consumes_ref = {elt.bam: elt.consumes_ref for elt in cigar_ops}

# Amount of reference sequence consumed by a cigar element
def ref_consumed(cigar_tuple):
    if len(cigar_tuple) != 2:
        raise ValueError("Invalid cigar tuple: %s" % cigar_tuple)
    if bam_to_consumes_ref[cigar_tuple[0]]:
        return cigar_tuple[1]
    else:
        return 0

# Total amount of reference sequence consumed by a list of cigar tuples
# Tuples are (bam, length)
def cigar_span(cigar_tuples):
    return sum(map(ref_consumed, cigar_tuples))

