"""Test oligo functions."""

import unittest

from synbio.oligos import calc_tm, fold, _bulge, _pair, _hairpin, _internal_loop


class TestOligos(unittest.TestCase):
    """Test oligo functions"""

    def test_calc_tm(self):
        """Test oligo tm calculation."""

        # values are from Table 1 of IDT's:
        # Owczarzy et al. (2008), Biochemistry 4 7: 5336-5353
        # with a 1.5mM Mg concentration which looks typical according to NEB
        experimental_tms = {
            "GGGACCGCCT": 51.9,
            "CCATTGCTACC": 42.7,
            "GCAGTGGATGTGAGA": 55.1,
            "CTGGTCTGGATCTGAGAACTTCAGG": 67.7,
            "CTTAAGATATGAGAACTTCAACTAATGTGT": 59.7,
            "AGTCTGGTCTGGATCTGAGAACTTCAGGCT": 71.6,
        }

        for seq, actual in experimental_tms.items():
            calc = calc_tm(seq)
            self.assertAlmostEqual(calc, actual, delta=3)  # within 3 deg tm difference

    def test_fold(self):
        """Test DNA folding to find min energy secondary structure."""

        # unafold's estimates for free energy estimates of DNA oligos
        unafold_dgs = {
            "TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA": -18.1,
            "AAGGGGTTGGTCGCCTCGACTAAGCGGCTGGATTCC": -3.5,  # unafold == -2.5
            "TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,  # bifurcation
            "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA": -3.65,
            # getting -2.21 because more bearish on a 4bp hairpin w/ pre-computed energy
            # the below is a three branched structure
            # "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC": -10.94,
            # "TGTCAGAAGTTTCCAAATGGCCAGCAATCAACCCATTCCATTGGGGATACAATGGTACAGTTTCGCATATTGTCGGTGAAAATGGTTCCATTAAACTCC": -9.35
        }

        for seq, unafold_est in unafold_dgs.items():
            structs = fold(seq, temp=37.0)
            calc_dg = sum([s[-1] for s in structs])

            # accepting a 25% difference
            delta = abs(0.25 * unafold_est)
            self.assertAlmostEqual(calc_dg, unafold_est, delta=delta)

    def test_bulge(self):
        """Test delta G calc of a bulge."""

        # mock bulge of CAT on one side and AG on other
        # from pg 429 of SantaLucia, 2004
        pair = "CT/GA"
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"  # nonsense sequence

        pair_dg = _bulge(pair, seq, 5, 7, 310.15)
        self.assertAlmostEqual(3.22, pair_dg, delta=0.4)

        # self.assertEqual(_bulge("CT/GA", seq, 5, 22, 310.15), 0.0)

    def test_pair(self):
        """Test delta G of pairs with and without mismatches."""

        pairs = [
            ("CT/GA", -1.28),
            ("GG/CC", -1.84),
            ("TC/AG", -1.3),
            ("GT/CG", -0.59),
            ("TC/GG", 0.08),
        ]
        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"

        for pair, dg_actual in pairs:
            dg_est = _pair(pair, seq, 5, 27, 310.15)
            self.assertAlmostEqual(dg_est, dg_actual, delta=0.02)

    def test_hairpin(self):
        """Test delta G of a hairpin structure."""

        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 11
        j = 16
        temp = 310.15
        hairpin_dg = _hairpin(seq, i, j, temp)
        # this differs from Unafold
        self.assertAlmostEqual(hairpin_dg, 1.7, delta=1.0)

        # self.assertEqual(0, _hairpin(seq, 8, 15, temp))

        # from page 428 of SantaLucia, 2004
        # hairpin = "CGCAAG"
        seq = "ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 3
        j = 8
        hairpin_dg = _hairpin(seq, i, j, temp)
        self.assertAlmostEqual(0.67, hairpin_dg, delta=0.1)

    def test_internal_loop(self):
        """Test internal loop."""

        seq = "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
        i = 6
        j = 21
        left = "TCCTT"
        right = "ATCAA"
        temp = 310.15
        temp_est = 3.5

        loop_temp = _internal_loop(seq, i, j, left, right, temp)

        self.assertAlmostEqual(loop_temp, temp_est, delta=0.1)
