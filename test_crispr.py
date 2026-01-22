import unittest
import crispr_designer as cd

class TestCRISPR(unittest.TestCase):

    def test_reverse_complement(self):
        self.assertEqual(cd.reverse_complement("ATCG"), "CGAT")

    def test_scoring(self):
        self.assertGreater(
            cd.score_grna("GCGCGCATATATGCGCGCAT"),
            cd.score_grna("AAAAAAAAAAAAAAAAAAAA")
        )

if __name__ == "__main__":
    unittest.main()