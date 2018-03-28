import unittest
from sam_utils import cigar_span

class Tests(unittest.TestCase):
    
    def test_cigar_len(self):
        self.assertEqual(cigar_span([(0, 10)]), 10)
        self.assertEqual(cigar_span([(1, 10)]), 0)
        self.assertEqual(cigar_span([(1, 100), (4, 100), (7, 10), (2, 10)]), 20)
        
        