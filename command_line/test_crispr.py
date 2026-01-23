import unittest
import crispr_designer as cd

class TestCRISPR(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEqual(cd.reverse_complement("ATCG"), "CGAT")

    def test_scoring(self):
        self.assertGreater(
            cd.score_grna("GCGCGCATATATGCGCGCAT", "GCGCGCATATATGCGCGCATCGATCGAT"),
            cd.score_grna("AAAAAAAAAAAAAAAAAAAA", "AAAAAAAAAAAAAAAAAAAAAAAAAAAA")
        )

    def test_hamming(self):
        assert cd.hamming_distance("AAAA", "AAAT") == 1

    def test_off_target_count(self):
        dna = "AAAAAAAAAAAAAAAAAAAA"
        guide = "AAAAAAAAAAAAAAAAAAAA"
        assert cd.count_off_targets(guide, dna) >= 0

if __name__ == "__main__":
    unittest.main()