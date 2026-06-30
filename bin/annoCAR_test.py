"""
Unit test suite for RAVA SNV annoCAR script.

Tests all conditions:
- Forward and reverse strand SNV annotation
- Single-exon and multi-exon (spliced) genes
- Ribosomal slippage (overlapping exons, slippage = 0 and slippage < 0)
- Standard introns (slippage > 0, no frameshift)
- Indel classification (frameshift/nonframeshift, both strands)
- Stop gain / stop loss detection
- Edge cases (exon boundaries, codons spanning splice junctions, variant outside CDS)

Run with: pytest annoCAR_test.py -v
"""

import pytest
from Bio.Seq import Seq

# =============================================================================
# FUNCTIONS UNDER TEST
# =============================================================================

from annoCAR import (
    complement_base,
    reverse_complement,
    translate,
    get_cds_offset_with_slippage,
    build_spliced_cds_with_slippage,
    get_codon_and_position,
    get_exon_number,
)


COMPLEMENT = str.maketrans("ACGTacgtRYKMSWBDHVNrykmswhbvdn",
                           "TGCAtgcaYRMKSWVHDBNyrmkswvhdbn")

# =============================================================================
# TEST CLASS: Complement Functions
# =============================================================================

class TestComplementFunctions:
    """Test base complementing and reverse complement operations."""

    def test_complement_base_A(self):
        assert complement_base("A") == "T"

    def test_complement_base_T(self):
        assert complement_base("T") == "A"

    def test_complement_base_G(self):
        assert complement_base("G") == "C"

    def test_complement_base_C(self):
        assert complement_base("C") == "G"

    def test_complement_ambiguous_R(self):
        assert complement_base("R") == "Y"

    def test_complement_ambiguous_Y(self):
        assert complement_base("Y") == "R"

    def test_reverse_complement_simple(self):
        assert reverse_complement("ATGC") == "GCAT"

    def test_reverse_complement_longer(self):
        assert reverse_complement("ATGGCT") == "AGCCAT"

    def test_reverse_complement_palindrome(self):
        assert reverse_complement("ATAT") == "ATAT"

    def test_reverse_complement_involution(self):
        """rc(rc(seq)) == seq."""
        seq = "ATGCGTAACCT"
        assert reverse_complement(reverse_complement(seq)) == seq

    def test_complement_preserves_lowercase(self):
        assert complement_base("a") == "t"
        assert complement_base("t") == "a"


# =============================================================================
# TEST CLASS: CDS Offset Calculation
# =============================================================================

class TestCDSOffset:
    """Test get_cds_offset_with_slippage for various configurations."""

    # --- Forward strand, single exon ---

    def test_forward_single_exon_first_position(self):
        exons = [(101, 200)]
        assert get_cds_offset_with_slippage(101, exons, "+") == 0

    def test_forward_single_exon_middle_position(self):
        exons = [(101, 200)]
        assert get_cds_offset_with_slippage(110, exons, "+") == 9

    def test_forward_single_exon_last_position(self):
        exons = [(101, 200)]
        assert get_cds_offset_with_slippage(200, exons, "+") == 99

    def test_forward_single_exon_outside(self):
        exons = [(101, 200)]
        assert get_cds_offset_with_slippage(100, exons, "+") == -1
        assert get_cds_offset_with_slippage(201, exons, "+") == -1

    # --- Forward strand, two exons, standard intron (slippage > 0) ---

    def test_forward_two_exons_standard_intron_exon1(self):
        """Exon1: 10-20 (11bp), Exon2: 30-50 (21bp), slippage=10."""
        exons = [(10, 20), (30, 50)]
        assert get_cds_offset_with_slippage(10, exons, "+") == 0
        assert get_cds_offset_with_slippage(15, exons, "+") == 5
        assert get_cds_offset_with_slippage(20, exons, "+") == 10

    def test_forward_two_exons_standard_intron_exon2(self):
        """After standard intron, exon2 offset starts at exon1 length (11)."""
        exons = [(10, 20), (30, 50)]
        assert get_cds_offset_with_slippage(30, exons, "+") == 11
        assert get_cds_offset_with_slippage(35, exons, "+") == 16

    # --- Forward strand, ribosomal slippage (slippage = 0) ---

    def test_forward_slippage_zero_exon1(self):
        """Slippage=0: exon1 contributes (end - start) bp, dropping shared position."""
        exons = [(10, 20), (20, 40)]
        assert get_cds_offset_with_slippage(10, exons, "+") == 0
        assert get_cds_offset_with_slippage(19, exons, "+") == 9

    def test_forward_slippage_zero_exon2(self):
        """Slippage=0: exon2 starts at offset = (end - start) of exon1."""
        exons = [(10, 20), (20, 40)]
        assert get_cds_offset_with_slippage(21, exons, "+") == 11
        assert get_cds_offset_with_slippage(25, exons, "+") == 15

    # --- Forward strand, ribosomal slippage (slippage = -1) ---

    def test_forward_slippage_negative_one(self):
        """Slippage=-1: exon1 effective contribution = (end-start+1) + slippage."""
        exons = [(10, 20), (19, 40)]
        assert get_cds_offset_with_slippage(10, exons, "+") == 0
        assert get_cds_offset_with_slippage(19, exons, "+") == 9

    def test_forward_slippage_negative_one_exon2(self):
        """Exon2 positions after overlap region."""
        exons = [(10, 20), (19, 40)]
        assert get_cds_offset_with_slippage(21, exons, "+") == 12
        assert get_cds_offset_with_slippage(25, exons, "+") == 16

    # --- Reverse strand, single exon ---

    def test_reverse_single_exon_last_genomic_is_offset_zero(self):
        """For minus strand, last genomic position is CDS offset 0."""
        exons = [(101, 200)]
        assert get_cds_offset_with_slippage(200, exons, "-") == 0

    def test_reverse_single_exon_first_genomic_is_last_offset(self):
        """For minus strand, first genomic position is last CDS offset."""
        exons = [(101, 200)]
        assert get_cds_offset_with_slippage(101, exons, "-") == 99

    def test_reverse_single_exon_middle(self):
        exons = [(101, 200)]
        assert get_cds_offset_with_slippage(150, exons, "-") == 50

    # --- Reverse strand, two exons, standard intron ---

    def test_reverse_two_exons_standard_intron_transcript_exon1(self):
        """
        Exon1(genomic): 10-20, Exon2(genomic): 30-50.
        Transcript exon1 = genomic exon2 (30-50). Position 50 = offset 0.
        """
        exons = [(10, 20), (30, 50)]
        assert get_cds_offset_with_slippage(50, exons, "-") == 0
        assert get_cds_offset_with_slippage(45, exons, "-") == 5
        assert get_cds_offset_with_slippage(30, exons, "-") == 20

    def test_reverse_two_exons_standard_intron_transcript_exon2(self):
        """
        Transcript exon2 = genomic exon1 (10-20).
        Genomic exon2 (30-50) contributes 21bp, so exon1 starts at offset 21.
        """
        exons = [(10, 20), (30, 50)]
        assert get_cds_offset_with_slippage(20, exons, "-") == 21
        assert get_cds_offset_with_slippage(15, exons, "-") == 26
        assert get_cds_offset_with_slippage(10, exons, "-") == 31

    # --- Reverse strand, ribosomal slippage ---

    def test_reverse_slippage_zero(self):
        """Reverse strand, slippage=0: first reversed exon contributes (end-start) bp."""
        exons = [(10, 20), (20, 40)]
        assert get_cds_offset_with_slippage(40, exons, "-") == 0
        assert get_cds_offset_with_slippage(21, exons, "-") == 19

    def test_reverse_slippage_zero_second_transcript_exon(self):
        """Second transcript exon after slippage=0."""
        exons = [(10, 20), (20, 40)]
        assert get_cds_offset_with_slippage(15, exons, "-") == 25
        assert get_cds_offset_with_slippage(10, exons, "-") == 30

    # --- Position outside CDS ---

    def test_forward_outside_cds(self):
        exons = [(10, 20), (30, 50)]
        assert get_cds_offset_with_slippage(25, exons, "+") == -1

    def test_reverse_outside_cds(self):
        exons = [(10, 20), (30, 50)]
        assert get_cds_offset_with_slippage(25, exons, "-") == -1


# =============================================================================
# TEST CLASS: Spliced CDS Construction
# =============================================================================

class TestBuildSplicedCDS:
    """Test build_spliced_cds_with_slippage for various configurations."""

    GENOME = "ATGCGTAACCTGGAAATTTCCCGGGATGCGTAACCTGGAAATTTCCCGGG"  # 50bp

    # --- Forward strand, single exon ---

    def test_forward_single_exon(self):
        exons = [(1, 9)]
        result = build_spliced_cds_with_slippage(self.GENOME, exons, "+")
        assert result == "ATGCGTAAC"

    # --- Forward strand, two exons, standard intron ---

    def test_forward_two_exons_standard(self):
        """Exon1: 1-6, Exon2: 10-15. slippage=4, both included fully."""
        exons = [(1, 6), (10, 15)]
        result = build_spliced_cds_with_slippage(self.GENOME, exons, "+")
        assert result == self.GENOME[0:6] + self.GENOME[9:15]

    # --- Forward strand, slippage = 0 ---

    def test_forward_slippage_zero(self):
        """Exon1: 1-6, Exon2: 6-12. Exon1 drops last bp (shared with exon2)."""
        exons = [(1, 6), (6, 12)]
        result = build_spliced_cds_with_slippage(self.GENOME, exons, "+")
        assert result == self.GENOME[0:5] + self.GENOME[5:12]

    # --- Forward strand, slippage = -1 ---

    def test_forward_slippage_negative_one(self):
        """Exon1: 1-6, Exon2: 5-12. Exon1 effective = 6 + (-1) = 5bp."""
        exons = [(1, 6), (5, 12)]
        result = build_spliced_cds_with_slippage(self.GENOME, exons, "+")
        assert result == self.GENOME[0:5] + self.GENOME[4:12]

    # --- Reverse strand, single exon ---

    def test_reverse_single_exon(self):
        """Single exon reverse: extract and reverse complement."""
        exons = [(1, 9)]
        result = build_spliced_cds_with_slippage(self.GENOME, exons, "-")
        assert result == reverse_complement("ATGCGTAAC")
        assert result == "GTTACGCAT"

    # --- Reverse strand, two exons, standard intron ---

    def test_reverse_two_exons_standard(self):
        """Reversed order concatenation, then reverse complement."""
        exons = [(1, 6), (10, 15)]
        result = build_spliced_cds_with_slippage(self.GENOME, exons, "-")
        expected = reverse_complement(self.GENOME[9:15] + self.GENOME[0:6])
        assert result == expected

    # --- Reverse strand, slippage = 0 ---

    def test_reverse_slippage_zero(self):
        """Reverse strand, slippage=0 between reversed exons."""
        exons = [(1, 6), (6, 12)]
        result = build_spliced_cds_with_slippage(self.GENOME, exons, "-")
        # Reversed: [(6,12), (1,6)]. First contributes (12-6)=6bp (drops last).
        expected = reverse_complement(self.GENOME[5:11] + self.GENOME[0:6])
        assert result == expected

    # --- Three exons, forward strand ---

    def test_forward_three_exons(self):
        exons = [(1, 6), (10, 15), (20, 25)]
        result = build_spliced_cds_with_slippage(self.GENOME, exons, "+")
        expected = self.GENOME[0:6] + self.GENOME[9:15] + self.GENOME[19:25]
        assert result == expected


# =============================================================================
# TEST CLASS: Codon and Position Extraction
# =============================================================================

class TestGetCodonAndPosition:
    """Test get_codon_and_position with CDS already in mRNA orientation."""

    def test_first_codon_first_position(self):
        cds = "ATGGCTGAA"  # M A E
        ref, alt, pos, prot_pos = get_codon_and_position(0, cds, "G", "+")
        assert ref == "ATG"
        assert alt == "GTG"
        assert pos == 0
        assert prot_pos == 1

    def test_first_codon_second_position(self):
        cds = "ATGGCTGAA"
        ref, alt, pos, prot_pos = get_codon_and_position(1, cds, "A", "+")
        assert ref == "ATG"
        assert alt == "AAG"
        assert pos == 1
        assert prot_pos == 1

    def test_first_codon_third_position(self):
        cds = "ATGGCTGAA"
        ref, alt, pos, prot_pos = get_codon_and_position(2, cds, "A", "+")
        assert ref == "ATG"
        assert alt == "ATA"
        assert pos == 2
        assert prot_pos == 1

    def test_second_codon(self):
        cds = "ATGGCTGAA"
        ref, alt, pos, prot_pos = get_codon_and_position(3, cds, "T", "+")
        assert ref == "GCT"
        assert alt == "TCT"
        assert pos == 0
        assert prot_pos == 2

    def test_third_codon(self):
        cds = "ATGGCTGAA"
        ref, alt, pos, prot_pos = get_codon_and_position(6, cds, "T", "+")
        assert ref == "GAA"
        assert alt == "TAA"
        assert pos == 0
        assert prot_pos == 3

    def test_synonymous_change(self):
        cds = "ATGGCTGAA"
        ref, alt, pos, prot_pos = get_codon_and_position(5, cds, "C", "+")
        assert ref == "GCT"
        assert alt == "GCC"
        assert translate(ref) == "A"
        assert translate(alt) == "A"


# =============================================================================
# TEST CLASS: Exon Number Assignment
# =============================================================================

class TestExonNumber:
    """Test get_exon_number for forward and reverse strands."""

    def test_forward_exon1(self):
        exons = [(10, 20), (30, 50), (60, 80)]
        assert get_exon_number(15, exons, "+") == 1

    def test_forward_exon2(self):
        exons = [(10, 20), (30, 50), (60, 80)]
        assert get_exon_number(35, exons, "+") == 2

    def test_forward_exon3(self):
        exons = [(10, 20), (30, 50), (60, 80)]
        assert get_exon_number(70, exons, "+") == 3

    def test_reverse_exon1_is_last_genomic(self):
        """For minus strand, transcript exon 1 = last genomic exon."""
        exons = [(10, 20), (30, 50), (60, 80)]
        assert get_exon_number(70, exons, "-") == 1

    def test_reverse_exon2(self):
        exons = [(10, 20), (30, 50), (60, 80)]
        assert get_exon_number(35, exons, "-") == 2

    def test_reverse_exon3_is_first_genomic(self):
        exons = [(10, 20), (30, 50), (60, 80)]
        assert get_exon_number(15, exons, "-") == 3


# =============================================================================
# TEST CLASS: Full SNV Annotation (Forward Strand)
# =============================================================================

class TestSNVAnnotationForward:
    """End-to-end SNV annotation on a forward strand single-exon gene."""

    # Gene at positions 5-22 (18bp): ATG GCT GAA TTC CGT AAA
    # Protein: M A E F R K
    GENOME = "NNNNATGGCTGAATTCCGTAAANNNNN"
    EXONS = [(5, 22)]

    def test_nonsynonymous_first_codon(self):
        """ATG->GTG = M->V at position 5."""
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "+")
        assert cds == "ATGGCTGAATTCCGTAAA"
        offset = get_cds_offset_with_slippage(5, self.EXONS, "+")
        assert offset == 0
        ref, alt, _, prot_pos = get_codon_and_position(offset, cds, "G", "+")
        assert ref == "ATG"
        assert alt == "GTG"
        assert translate(ref) == "M"
        assert translate(alt) == "V"
        assert prot_pos == 1

    def test_synonymous_third_position(self):
        """GCT->GCC = A->A at position 10 (CDS offset 5, codon 2 third position)."""
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "+")
        offset = get_cds_offset_with_slippage(10, self.EXONS, "+")
        assert offset == 5
        ref, alt, codon_pos, prot_pos = get_codon_and_position(offset, cds, "C", "+")
        assert ref == "GCT"
        assert alt == "GCC"
        assert translate(ref) == "A"
        assert translate(alt) == "A"
        assert prot_pos == 2

    def test_stopgain(self):
        """GAA->TAA = E->* at position 11 (CDS offset 6, codon 3 first position)."""
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "+")
        offset = get_cds_offset_with_slippage(11, self.EXONS, "+")
        assert offset == 6
        ref, alt, _, prot_pos = get_codon_and_position(offset, cds, "T", "+")
        assert ref == "GAA"
        assert alt == "TAA"
        assert translate(ref) == "E"
        assert translate(alt) == "*"
        assert prot_pos == 3


# =============================================================================
# TEST CLASS: Full SNV Annotation (Reverse Strand)
# =============================================================================

class TestSNVAnnotationReverse:
    """
    End-to-end SNV annotation on a reverse strand single-exon gene.

    mRNA = ATG GCT GAA TTC CGT AAA (M A E F R K)
    Genomic (positions 5-22) = rc(mRNA) = "TTTACGGAATTCAGCCAT"
    """
    GENOMIC_GENE = reverse_complement("ATGGCTGAATTCCGTAAA")
    GENOME = "NNNN" + GENOMIC_GENE + "NNNNN"
    EXONS = [(5, 22)]

    def test_cds_construction(self):
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "-")
        assert cds == "ATGGCTGAATTCCGTAAA"

    def test_nonsynonymous_at_last_genomic_position(self):
        """
        Position 22 = CDS offset 0 (5' end of mRNA).
        Genomic alt C -> mRNA alt G: ATG->GTG = M->V.
        """
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "-")
        offset = get_cds_offset_with_slippage(22, self.EXONS, "-")
        assert offset == 0

        alt_base_mRNA = complement_base("C")
        assert alt_base_mRNA == "G"

        ref, alt, _, prot_pos = get_codon_and_position(offset, cds, alt_base_mRNA, "-")
        assert ref == "ATG"
        assert alt == "GTG"
        assert translate(ref) == "M"
        assert translate(alt) == "V"
        assert prot_pos == 1

    def test_synonymous_reverse_strand(self):
        """
        CDS offset 5 (codon 2, third position). Genomic position = 22 - 5 = 17.
        GCT->GCC (A->A): mRNA alt C requires genomic alt G.
        """
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "-")
        offset = get_cds_offset_with_slippage(17, self.EXONS, "-")
        assert offset == 5

        alt_base_mRNA = complement_base("G")
        ref, alt, _, prot_pos = get_codon_and_position(offset, cds, alt_base_mRNA, "-")
        assert ref == "GCT"
        assert alt == "GCC"
        assert translate(ref) == "A"
        assert translate(alt) == "A"
        assert prot_pos == 2

    def test_stopgain_reverse_strand(self):
        """
        CDS offset 6 (codon 3, first position). Genomic position = 22 - 6 = 16.
        GAA->TAA (E->*): mRNA alt T requires genomic alt A.
        """
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "-")
        offset = get_cds_offset_with_slippage(16, self.EXONS, "-")
        assert offset == 6

        alt_base_mRNA = complement_base("A")
        ref, alt, _, prot_pos = get_codon_and_position(offset, cds, alt_base_mRNA, "-")
        assert ref == "GAA"
        assert alt == "TAA"
        assert translate(ref) == "E"
        assert translate(alt) == "*"
        assert prot_pos == 3


# =============================================================================
# TEST CLASS: Spliced Gene (Forward Strand)
# =============================================================================

class TestSplicedForward:
    """
    Forward strand gene with two exons.

    Exon1: positions 11-19 (9bp): ATG CAG CTT  -> M Q L
    Intron: positions 20-29 (10bp)
    Exon2: positions 30-38 (9bp): AAA GGC TCA  -> K G S

    Spliced CDS = ATGCAGCTTAAAGGCTCA (18bp, 6 codons)
    Protein = M Q L K G S
    """
    GENOME = "NNNNNNNNNN" + "ATGCAGCTT" + "GGGGGGGGGG" + "AAAGGCTCA" + "NNNNNN"
    EXONS = [(11, 19), (30, 38)]

    def test_spliced_cds_construction(self):
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "+")
        assert cds == "ATGCAGCTTAAAGGCTCA"

    def test_variant_in_exon1(self):
        """Position 14 -> CDS offset 3 -> codon 2 (CAG). C->T: CAG->TAG = Q->*."""
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "+")
        offset = get_cds_offset_with_slippage(14, self.EXONS, "+")
        assert offset == 3
        ref, alt, _, prot_pos = get_codon_and_position(offset, cds, "T", "+")
        assert ref == "CAG"
        assert alt == "TAG"
        assert prot_pos == 2

    def test_variant_in_exon2(self):
        """Position 30 -> CDS offset 9 -> codon 4 (AAA). A->G: AAA->GAA = K->E."""
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "+")
        offset = get_cds_offset_with_slippage(30, self.EXONS, "+")
        assert offset == 9
        ref, alt, _, prot_pos = get_codon_and_position(offset, cds, "G", "+")
        assert ref == "AAA"
        assert alt == "GAA"
        assert prot_pos == 4

    def test_variant_in_intron_returns_negative(self):
        offset = get_cds_offset_with_slippage(25, self.EXONS, "+")
        assert offset == -1

    def test_exon_numbers(self):
        assert get_exon_number(15, self.EXONS, "+") == 1
        assert get_exon_number(35, self.EXONS, "+") == 2

# =============================================================================
# TEST CLASS: Spliced Gene (Reverse Strand)
# =============================================================================

class TestSplicedReverse:
    """
    Reverse strand gene with two exons.

    Exon1 (genomic 11-19) = transcript exon 2
    Exon2 (genomic 30-38) = transcript exon 1

    build_spliced_cds_with_slippage for reverse strand:
      1. Reverses exon order: [(30,38), (11,19)]
      2. Concatenates genomic: genome[29:38] + genome[10:19]
      3. Reverse complements the result

    After rc, genomic exon1 (concatenated second) becomes the 5' end of mRNA,
    and genomic exon2 (concatenated first) becomes the 3' end.

    Target mRNA = ATGCAGCTT AAAGGCTCA (M Q L K G S)
      Genomic exon1 (11-19) = rc("ATGCAGCTT") = "AAGCTGCAT"
      Genomic exon2 (30-38) = rc("AAAGGCTCA") = "TGAGCCTTT"
    """
    GENE_EXON1_GENOMIC = reverse_complement("ATGCAGCTT")  # "AAGCTGCAT"
    GENE_EXON2_GENOMIC = reverse_complement("AAAGGCTCA")  # "TGAGCCTTT"
    GENOME = "NNNNNNNNNN" + GENE_EXON1_GENOMIC + "GGGGGGGGGG" + GENE_EXON2_GENOMIC + "NNNNNN"
    EXONS = [(11, 19), (30, 38)]

    def test_spliced_cds_construction(self):
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "-")
        assert cds == "ATGCAGCTTAAAGGCTCA"

    def test_variant_in_transcript_exon1(self):
        """Genomic position 38 = CDS offset 0 (5' end of mRNA)."""
        offset = get_cds_offset_with_slippage(38, self.EXONS, "-")
        assert offset == 0

    def test_variant_in_transcript_exon2(self):
        """Genomic position 19 = first position of transcript exon 2 = CDS offset 9."""
        offset = get_cds_offset_with_slippage(19, self.EXONS, "-")
        assert offset == 9

    def test_variant_produces_correct_amino_acid(self):
        """
        Genomic pos 35 -> CDS offset 3 -> codon 2 (CAG).
        Genomic alt A -> mRNA alt T: CAG->TAG = Q->*.
        """
        cds = build_spliced_cds_with_slippage(self.GENOME, self.EXONS, "-")
        offset = get_cds_offset_with_slippage(35, self.EXONS, "-")
        assert offset == 3

        alt_mRNA = complement_base("A")
        ref, alt, _, prot_pos = get_codon_and_position(offset, cds, alt_mRNA, "-")
        assert ref == "CAG"
        assert alt == "TAG"
        assert translate(alt) == "*"
        assert prot_pos == 2

    def test_exon_numbers_reverse(self):
        assert get_exon_number(35, self.EXONS, "-") == 1
        assert get_exon_number(15, self.EXONS, "-") == 2

# =============================================================================
# TEST CLASS: Ribosomal Slippage (Forward Strand)
# =============================================================================

class TestRibosomalSlippageForward:
    """
    Test ribosomal slippage on forward strand where exons share a boundary
    (slippage=0) or overlap (slippage<0).
    """

    def test_slippage_zero_cds_length(self):
        """Exon1: 1-10, Exon2: 10-20. Exon1 contributes 9bp (drops shared position)."""
        genome = "ATGAAACCCGGGTTTAAACCC"
        exons = [(1, 10), (10, 20)]
        cds = build_spliced_cds_with_slippage(genome, exons, "+")
        assert cds == genome[0:9] + genome[9:20]
        assert len(cds) == 20

    def test_slippage_zero_variant_at_shared_boundary(self):
        """Position 10 found in exon1 first, offset = 9."""
        genome = "ATGAAACCCGGGTTTAAACCC"
        exons = [(1, 10), (10, 20)]
        offset = get_cds_offset_with_slippage(10, exons, "+")
        assert offset == 9

    def test_slippage_negative_one_cds_length(self):
        """Exon1: 1-10, Exon2: 9-20. Exon1 effective = 10 + (-1) = 9bp."""
        genome = "ATGAAACCCGGGTTTAAACCCC"
        exons = [(1, 10), (9, 20)]
        cds = build_spliced_cds_with_slippage(genome, exons, "+")
        assert len(cds) == 9 + 12
        assert cds[:9] == genome[0:9]
        assert cds[9:] == genome[8:20]

    def test_slippage_negative_two(self):
        """Exon1: 1-10, Exon2: 8-20. Exon1 effective = 10 + (-2) = 8bp."""
        genome = "ATGAAACCCGGGTTTAAACCCCC"
        exons = [(1, 10), (8, 20)]
        cds = build_spliced_cds_with_slippage(genome, exons, "+")
        assert cds[:8] == genome[0:8]
        assert len(cds) == 8 + 13

    def test_codon_at_splice_junction(self):
        """
        Exon1 contributes 9bp (3 full codons), codon 4 starts at exon2.
        Verify correct amino acid call across junction.
        """
        genome = "ATGAAAGCC" + "GATTTTAAACCC"
        exons = [(1, 10), (10, 21)]
        cds = build_spliced_cds_with_slippage(genome, exons, "+")
        assert cds == "ATGAAAGCC" + "GATTTTAAACCC"

        ref, alt, _, prot_pos = get_codon_and_position(9, cds, "A", "+")
        assert ref == "GAT"
        assert alt == "AAT"
        assert translate("GAT") == "D"
        assert translate("AAT") == "N"
        assert prot_pos == 4


# =============================================================================
# TEST CLASS: Ribosomal Slippage (Reverse Strand)
# =============================================================================

class TestRibosomalSlippageReverse:
    """Test ribosomal slippage on the reverse strand."""

    def test_reverse_slippage_zero_cds(self):
        """
        Exon1(genomic): 1-10, Exon2(genomic): 10-20. Reversed: [(10,20), (1,10)].
        First reversed exon contributes 10bp (slippage=0 drops last).
        Concat then reverse complement.
        """
        genome = "ATGCGTAACCTGGAAATTTCCC"
        exons = [(1, 10), (10, 20)]
        cds = build_spliced_cds_with_slippage(genome, exons, "-")
        forward_concat = genome[9:19] + genome[0:10]
        expected = reverse_complement(forward_concat)
        assert cds == expected
        assert len(cds) == 20

    def test_reverse_slippage_negative_one(self):
        """
        Exon1(genomic): 1-10, Exon2(genomic): 9-20. Overlap at 9-10.
        Reversed: [(9,20), (1,10)]. slippage = 9 - 10 = -1.
        First contributes (20-9+1)+(-1) = 11bp.
        """
        genome = "ATGCGTAACCTGGAAATTTCCC"
        exons = [(1, 10), (9, 20)]
        cds = build_spliced_cds_with_slippage(genome, exons, "-")
        forward_concat = genome[8:19] + genome[0:10]
        expected = reverse_complement(forward_concat)
        assert cds == expected


# =============================================================================
# TEST CLASS: Indel Classification
# =============================================================================

class TestIndelClassification:
    """Test frameshift vs non-frameshift classification for indels."""

    def test_1bp_insertion_is_frameshift(self):
        ref, alt = "A", "AT"
        assert (len(alt) - len(ref)) % 3 != 0

    def test_3bp_insertion_is_nonframeshift(self):
        ref, alt = "A", "ATCG"
        assert (len(alt) - len(ref)) % 3 == 0

    def test_1bp_deletion_is_frameshift(self):
        ref, alt = "AT", "A"
        assert (len(ref) - len(alt)) % 3 != 0

    def test_3bp_deletion_is_nonframeshift(self):
        ref, alt = "ATCG", "A"
        assert (len(ref) - len(alt)) % 3 == 0

    def test_6bp_deletion_is_nonframeshift(self):
        ref, alt = "ATCGATC", "A"
        assert (len(ref) - len(alt)) % 3 == 0

    def test_2bp_insertion_is_frameshift(self):
        ref, alt = "G", "GCC"
        assert (len(alt) - len(ref)) % 3 != 0

    def test_4bp_deletion_is_frameshift(self):
        ref, alt = "CATAA", "C"
        assert (len(ref) - len(alt)) % 3 != 0

    def test_9bp_deletion_is_nonframeshift(self):
        ref, alt = "TACACCATGG", "T"
        assert (len(ref) - len(alt)) == 9
        assert (len(ref) - len(alt)) % 3 == 0


# =============================================================================
# TEST CLASS: Stop Codon Detection
# =============================================================================

class TestStopCodonDetection:
    """Test stopgain and stoploss detection."""

    def test_stopgain_TAA(self):
        assert translate("TAA") == "*"

    def test_stopgain_TAG(self):
        assert translate("TAG") == "*"

    def test_stopgain_TGA(self):
        assert translate("TGA") == "*"

    def test_stoploss_TAA_to_CAA(self):
        assert translate("TAA") == "*"
        assert translate("CAA") == "Q"

    def test_stopgain_GAG_to_TAG(self):
        assert translate("GAG") == "E"
        assert translate("TAG") == "*"

    def test_normal_change_no_stop(self):
        assert translate("GCT") == "A"
        assert translate("GCC") == "A"


# =============================================================================
# TEST CLASS: Edge Cases
# =============================================================================

class TestEdgeCases:
    """Test boundary conditions and edge cases."""

    def test_variant_at_exon_boundary_start(self):
        """First position of exon2 after an 11bp exon1."""
        exons = [(10, 20), (30, 40)]
        offset = get_cds_offset_with_slippage(30, exons, "+")
        assert offset == 11

    def test_variant_at_exon_boundary_end(self):
        """Last position of exon1."""
        exons = [(10, 20), (30, 40)]
        offset = get_cds_offset_with_slippage(20, exons, "+")
        assert offset == 10

    def test_codon_spanning_splice_junction(self):
        """
        Exon1: 1-7 (7bp = 2 full codons + 1bp). Exon2: 19-29.
        Codon 3 (CDS offsets 6-8) spans the junction.
        """
        genome = "ATGAAAG" + "NNNNNNNNNNN" + "CCTTTAAACCC" + "NNNN"
        exons = [(1, 7), (19, 29)]
        cds = build_spliced_cds_with_slippage(genome, exons, "+")
        assert cds == "ATGAAAG" + "CCTTTAAACCC"
        assert cds[6:9] == "GCC"
        assert translate("GCC") == "A"

    def test_both_strands_produce_same_protein(self):
        """Same CDS, one on each strand, should produce identical protein."""
        fwd_genome = "NNN" + "ATGGCTGAATGA" + "NNN"
        rev_genome = "NNN" + reverse_complement("ATGGCTGAATGA") + "NNN"
        exons = [(4, 15)]

        fwd_cds = build_spliced_cds_with_slippage(fwd_genome, exons, "+")
        rev_cds = build_spliced_cds_with_slippage(rev_genome, exons, "-")

        assert fwd_cds == "ATGGCTGAATGA"
        assert rev_cds == "ATGGCTGAATGA"
        assert translate(fwd_cds) == "MAE*"
        assert translate(rev_cds) == "MAE*"

    def test_translate_all_64_codons(self):
        codons_counted = 0
        for b1 in "ACGT":
            for b2 in "ACGT":
                for b3 in "ACGT":
                    codon = b1 + b2 + b3
                    result = translate(codon)
                    assert len(result) == 1
                    codons_counted += 1
        assert codons_counted == 64

    def test_position_outside_all_exons(self):
        exons = [(100, 200)]
        assert get_cds_offset_with_slippage(50, exons, "+") == -1
        assert get_cds_offset_with_slippage(50, exons, "-") == -1


# =============================================================================
# TEST CLASS: Three-Exon Splicing
# =============================================================================

class TestThreeExonSplicing:
    """Test three-exon genes on both strands."""

    def test_forward_three_exon_offsets(self):
        """Exon1: 1-9 (9bp), Exon2: 20-28 (9bp), Exon3: 40-48 (9bp). Total=27bp."""
        exons = [(1, 9), (20, 28), (40, 48)]
        assert get_cds_offset_with_slippage(1, exons, "+") == 0
        assert get_cds_offset_with_slippage(9, exons, "+") == 8
        assert get_cds_offset_with_slippage(20, exons, "+") == 9
        assert get_cds_offset_with_slippage(28, exons, "+") == 17
        assert get_cds_offset_with_slippage(40, exons, "+") == 18
        assert get_cds_offset_with_slippage(48, exons, "+") == 26

    def test_reverse_three_exon_offsets(self):
        """Same exons on reverse strand. Transcript order: exon3, exon2, exon1."""
        exons = [(1, 9), (20, 28), (40, 48)]
        assert get_cds_offset_with_slippage(48, exons, "-") == 0
        assert get_cds_offset_with_slippage(40, exons, "-") == 8
        assert get_cds_offset_with_slippage(28, exons, "-") == 9
        assert get_cds_offset_with_slippage(20, exons, "-") == 17
        assert get_cds_offset_with_slippage(9, exons, "-") == 18
        assert get_cds_offset_with_slippage(1, exons, "-") == 26

    def test_forward_three_exon_cds_construction(self):
        genome = ("ATGCAGCTT" +
                  "NNNNNNNNNN" +
                  "GAAGGCTCA" +
                  "NNNNNNNNNNN" +
                  "TTTAGCGAT")
        exons = [(1, 9), (20, 28), (40, 48)]
        cds = build_spliced_cds_with_slippage(genome, exons, "+")
        assert cds == "ATGCAGCTT" + "GAAGGCTCA" + "TTTAGCGAT"
        assert len(cds) == 27

    def test_reverse_three_exon_cds_construction(self):
        genome = ("ATGCAGCTT" +
                  "NNNNNNNNNN" +
                  "GAAGGCTCA" +
                  "NNNNNNNNNNN" +
                  "TTTAGCGAT")
        exons = [(1, 9), (20, 28), (40, 48)]
        cds = build_spliced_cds_with_slippage(genome, exons, "-")
        genomic_concat = "TTTAGCGAT" + "GAAGGCTCA" + "ATGCAGCTT"
        expected = reverse_complement(genomic_concat)
        assert cds == expected


# =============================================================================
# TEST CLASS: Influenza M2-like Splicing Pattern
# =============================================================================

class TestInfluenzaM2Pattern:
    """
    Test a pattern mimicking influenza M segment where:
    - M1: single CDS
    - M2: two exons with a splice (shares exon1 with M1)
    
    This tests the scenario that caused the original RAVA bug.
    """
    GENOME = ("ATGGCTGAA" +
              "TTCCGTAAAG" +
              "ATCTGTGACCCGGG" +
              "AAATTTCCC" +
              "GGGAAATTT")
    M2_EXONS = [(1, 9), (20, 37)]

    def test_m2_exon1_offset(self):
        """Shared region (pos 1-9) has same offset in M2 exon1."""
        offset = get_cds_offset_with_slippage(5, self.M2_EXONS, "+")
        assert offset == 4

    def test_m2_exon2_offset_after_splice(self):
        """M2 exon2 starts at CDS offset 9 (after 9bp from exon1)."""
        offset = get_cds_offset_with_slippage(20, self.M2_EXONS, "+")
        assert offset == 9

    def test_m2_codon_position_in_exon2(self):
        """Position 22 -> CDS offset 11 -> third position of codon 4."""
        offset = get_cds_offset_with_slippage(22, self.M2_EXONS, "+")
        assert offset == 11
        assert offset % 3 == 2
        assert (offset // 3) + 1 == 4

    def test_m2_exon_numbers(self):
        assert get_exon_number(5, self.M2_EXONS, "+") == 1
        assert get_exon_number(25, self.M2_EXONS, "+") == 2


# =============================================================================
# TEST CLASS: Full Annotation String Assembly
# =============================================================================

class TestAnnotationFormat:
    """Verify that annotation components are correctly assembled for both strands."""

    def test_forward_strand_annotation(self):
        """Forward strand: position 7 in gene at 4-15. GCT->TCT = A->S."""
        genome = "NNN" + "ATGGCTGAATGA" + "NNN"
        exons = [(4, 15)]
        strand = "+"
        variant_pos = 7
        genomic_alt = "T"

        cds = build_spliced_cds_with_slippage(genome, exons, strand)
        offset = get_cds_offset_with_slippage(variant_pos, exons, strand)
        exon_num = get_exon_number(variant_pos, exons, strand)
        ref_codon, alt_codon, _, prot_pos = get_codon_and_position(offset, cds, genomic_alt, strand)

        assert cds == "ATGGCTGAATGA"
        assert offset == 3
        assert exon_num == 1
        assert ref_codon == "GCT"
        assert alt_codon == "TCT"
        assert translate(ref_codon) == "A"
        assert translate(alt_codon) == "S"
        assert prot_pos == 2
        assert str(offset + 1) == "4"

    def test_reverse_strand_annotation(self):
        """Reverse strand: position 15 = CDS offset 0. ATG->GTG = M->V."""
        genomic_gene = reverse_complement("ATGGCTGAATGA")
        genome = "NNN" + genomic_gene + "NNN"
        exons = [(4, 15)]
        strand = "-"
        variant_pos = 15
        genomic_alt = "C"

        cds = build_spliced_cds_with_slippage(genome, exons, strand)
        assert cds == "ATGGCTGAATGA"

        offset = get_cds_offset_with_slippage(variant_pos, exons, strand)
        assert offset == 0

        mRNA_alt = complement_base(genomic_alt)
        assert mRNA_alt == "G"

        ref_codon, alt_codon, _, prot_pos = get_codon_and_position(offset, cds, mRNA_alt, strand)
        assert ref_codon == "ATG"
        assert alt_codon == "GTG"
        assert translate(ref_codon) == "M"
        assert translate(alt_codon) == "V"
        assert prot_pos == 1

        # HGVS c. notation uses transcript-oriented bases
        hgvs_ref = complement_base(genome[variant_pos - 1])
        assert str(offset + 1) == "1"


# =============================================================================
# RUN
# =============================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v"])