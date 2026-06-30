import pandas as pd
import io
import os
from Bio import SeqIO
from Bio.Seq import Seq
import re
import warnings
import argparse
import sys


warnings.simplefilter(action="ignore", category=FutureWarning)


# =============================================================================
# COMPLEMENT TABLE
# =============================================================================

COMPLEMENT = str.maketrans("ACGTacgtRYKMSWBDHVNrykmswhbvdn",
                           "TGCAtgcaYRMKSWVHDBNyrmkswvhdbn")


def complement_base(base):
    """Complement a single nucleotide base."""
    return base.translate(COMPLEMENT)


def reverse_complement(seq):
    """Reverse complement a nucleotide sequence string."""
    return seq.translate(COMPLEMENT)[::-1]


# =============================================================================
# UTILITY FUNCTIONS
# =============================================================================

def read_vcf(path):
    with open(path, "r") as f:
        lines = [l for l in f if not l.startswith("##")]
    return pd.read_csv(
        io.StringIO("".join(lines)),
        dtype={
            "#CHROM": str,
            "POS": int,
            "ID": str,
            "REF": str,
            "ALT": str,
            "QUAL": str,
            "FILTER": str,
            "INFO": str,
        },
        sep="\t",
    ).rename(columns={"#CHROM": "CHROM"})


def translate(seq):
    table = {
        "ATA": "I", "ATC": "I", "ATT": "I", "ATG": "M",
        "ACA": "T", "ACC": "T", "ACG": "T", "ACT": "T",
        "AAC": "N", "AAT": "N", "AAA": "K", "AAG": "K",
        "AGC": "S", "AGT": "S", "AGA": "R", "AGG": "R",
        "CTA": "L", "CTC": "L", "CTG": "L", "CTT": "L",
        "CCA": "P", "CCC": "P", "CCG": "P", "CCT": "P",
        "CAC": "H", "CAT": "H", "CAA": "Q", "CAG": "Q",
        "CGA": "R", "CGC": "R", "CGG": "R", "CGT": "R",
        "GTA": "V", "GTC": "V", "GTG": "V", "GTT": "V",
        "GCA": "A", "GCC": "A", "GCG": "A", "GCT": "A",
        "GAC": "D", "GAT": "D", "GAA": "E", "GAG": "E",
        "GGA": "G", "GGC": "G", "GGG": "G", "GGT": "G",
        "TCA": "S", "TCC": "S", "TCG": "S", "TCT": "S",
        "TTC": "F", "TTT": "F", "TTA": "L", "TTG": "L",
        "TAC": "Y", "TAT": "Y", "TAA": "*", "TAG": "*",
        "TGC": "C", "TGT": "C", "TGA": "*", "TGG": "W",
    }
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(seq), 3):
            codon = seq[i : i + 3]
            protein += table[codon]
    return protein


def degenerative(codon):
    standard_bases = set("ACGT")
    if all(c in standard_bases for c in codon):
        return codon

    if "R" in codon:
        return codonReplace_recursive(codon, "R", ["G", "A"])
    elif "Y" in codon:
        return codonReplace_recursive(codon, "Y", ["C", "T"])
    elif "K" in codon:
        return codonReplace_recursive(codon, "K", ["G", "T"])
    elif "M" in codon:
        return codonReplace_recursive(codon, "M", ["A", "C"])
    elif "W" in codon:
        return codonReplace_recursive(codon, "W", ["A", "T"])
    elif "S" in codon:
        return codonReplace_recursive(codon, "S", ["G", "C"])
    elif "B" in codon:
        return codonReplace_recursive(codon, "B", ["G", "T", "C"])
    elif "D" in codon:
        return codonReplace_recursive(codon, "D", ["G", "A", "T"])
    elif "H" in codon:
        return codonReplace_recursive(codon, "H", ["A", "C", "T"])
    elif "V" in codon:
        return codonReplace_recursive(codon, "V", ["G", "C", "A"])
    elif "N" in codon:
        return codonReplace_recursive(codon, "N", ["A", "G", "C", "T"])

    return codon


def codonReplace_recursive(codon, degen, options):
    translations = set()
    for base in options:
        resolved = codon.replace(degen, base, 1)
        resolved = degenerative(resolved)
        if len(resolved) == 3 and all(c in "ACGT" for c in resolved):
            resolved = translate(resolved)
            translations.add(resolved)

    if len(translations) == 1:
        return translations.pop()
    else:
        return "or".join(sorted(translations))


def extract_gene_name(attr_string):
    attr = str(attr_string).strip()
    for field in attr.split(";"):
        field = field.strip()
        if field.startswith("ID=gene:"):
            return field.replace("ID=gene:", "")
        elif field.startswith("ID=CDS:"):
            return field.replace("ID=CDS:", "")
    if "ID=" in attr:
        id_part = attr.split("ID=")[1].split(";")[0]
        if ":" in id_part:
            return id_part.split(":")[-1]
        return id_part
    return attr


# =============================================================================
# GENE/CDS MAP BUILDING
# =============================================================================

def build_gene_cds_map(gff):
    gene_cds = {}
    gene_strand = {}
    gene_coords = {}

    for i in range(len(gff.index)):
        if pd.isnull(gff.iloc[i, 2]):
            continue
        if gff.iloc[i, 2] == "CDS":
            cds_start = int(gff.iloc[i, 3])
            cds_end = int(gff.iloc[i, 4])
            cds_attr = str(gff.iloc[i, 8]).strip()
            gene_name = extract_gene_name(cds_attr)
            cds_strand = gff.iloc[i, 6]
            if gene_name not in gene_cds:
                gene_cds[gene_name] = []
                gene_strand[gene_name] = cds_strand
            gene_cds[gene_name].append((cds_start, cds_end))

    for i in range(len(gff.index)):
        if pd.isnull(gff.iloc[i, 2]):
            continue
        if gff.iloc[i, 2] == "gene":
            gene_attr = str(gff.iloc[i, 8]).strip()
            gene_name = extract_gene_name(gene_attr)
            if gene_name not in gene_strand:
                gene_strand[gene_name] = gff.iloc[i, 6]
            start = int(gff.iloc[i, 3])
            end = int(gff.iloc[i, 4])
            if gene_name in gene_coords:
                old_start, old_end = gene_coords[gene_name]
                gene_coords[gene_name] = (min(old_start, start), max(old_end, end))
            else:
                gene_coords[gene_name] = (start, end)

    for gene_name in gene_cds:
        gene_cds[gene_name].sort(key=lambda x: x[0])

    return gene_cds, gene_strand, gene_coords


# =============================================================================
# CDS OFFSET CALCULATION - STRAND AWARE
# =============================================================================

def get_cds_offset_with_slippage(variant_pos, cds_exons, strand):
    """
    Calculate the 0-based offset within the spliced CDS for a given genomic position.
    
    For forward strand (+): exons are processed 5'->3' (low to high coordinates).
    For reverse strand (-): exons are processed 3'->5' (high to low coordinates).
    
    Returns -1 if variant_pos is not within any CDS exon.
    """
    # Determine which exon contains the variant
    containing_idx = -1
    for idx, (exon_start, exon_end) in enumerate(cds_exons):
        if variant_pos >= exon_start and variant_pos <= exon_end:
            containing_idx = idx
            break

    if containing_idx == -1:
        return -1

    if strand == "+":
        # Forward strand: process exons in genomic order (already sorted low->high)
        cumulative = 0
        for idx in range(containing_idx):
            exon_start, exon_end = cds_exons[idx]
            next_start = cds_exons[idx + 1][0]
            slippage = next_start - exon_end

            if slippage >= 1:
                cumulative += (exon_end - exon_start + 1)
            elif slippage == 0:
                cumulative += (exon_end - exon_start)
            else:
                cumulative += (exon_end - exon_start + 1 + slippage)

        offset_in_exon = variant_pos - cds_exons[containing_idx][0]
        return cumulative + offset_in_exon

    else:
        # Reverse strand: CDS is read from the last exon backward
        # Reverse the exon order for CDS construction
        reversed_exons = list(reversed(cds_exons))
        # Find which index in reversed order contains our variant
        rev_containing_idx = -1
        for idx, (exon_start, exon_end) in enumerate(reversed_exons):
            if variant_pos >= exon_start and variant_pos <= exon_end:
                rev_containing_idx = idx
                break

        if rev_containing_idx == -1:
            return -1

        cumulative = 0
        for idx in range(rev_containing_idx):
            exon_start, exon_end = reversed_exons[idx]
            # For reverse strand, slippage between consecutive reversed exons
            if idx + 1 < len(reversed_exons):
                prev_exon_start, prev_exon_end = reversed_exons[idx + 1]
                slippage = exon_start - prev_exon_end
            else:
                slippage = 1  # default no overlap

            if slippage >= 1:
                cumulative += (exon_end - exon_start + 1)
            elif slippage == 0:
                cumulative += (exon_end - exon_start)
            else:
                cumulative += (exon_end - exon_start + 1 + slippage)

        # For reverse strand, offset within the exon is from the END
        exon_start, exon_end = reversed_exons[rev_containing_idx]
        offset_in_exon = exon_end - variant_pos
        return cumulative + offset_in_exon


# =============================================================================
# SPLICED CDS CONSTRUCTION - STRAND AWARE
# =============================================================================

def build_spliced_cds_with_slippage(fasta_seq, cds_exons, strand):
    """
    Build the spliced CDS nucleotide sequence from genomic FASTA.
    
    For forward strand: concatenate exons in genomic order.
    For reverse strand: concatenate exons in reverse genomic order, then
                        reverse complement the entire result.
    
    The returned sequence is always in mRNA orientation (5'->3' of transcript).
    """
    if strand == "+":
        spliced = ""
        for idx, (exon_start, exon_end) in enumerate(cds_exons):
            if idx == len(cds_exons) - 1:
                spliced += fasta_seq[exon_start - 1 : exon_end]
            else:
                next_start = cds_exons[idx + 1][0]
                slippage = next_start - exon_end

                if slippage >= 1:
                    spliced += fasta_seq[exon_start - 1 : exon_end]
                elif slippage == 0:
                    spliced += fasta_seq[exon_start - 1 : exon_end - 1]
                else:
                    effective_length = (exon_end - exon_start + 1) + slippage
                    spliced += fasta_seq[exon_start - 1 : exon_start - 1 + effective_length]
        return spliced

    else:
        # Reverse strand: process exons in reverse genomic order
        reversed_exons = list(reversed(cds_exons))
        spliced = ""
        for idx, (exon_start, exon_end) in enumerate(reversed_exons):
            if idx == len(reversed_exons) - 1:
                spliced += fasta_seq[exon_start - 1 : exon_end]
            else:
                # For reverse strand, slippage is between current and next
                # (next in reversed order = previous in genomic order)
                next_exon_start, next_exon_end = reversed_exons[idx + 1]
                slippage = exon_start - next_exon_end

                if slippage >= 1:
                    spliced += fasta_seq[exon_start - 1 : exon_end]
                elif slippage == 0:
                    spliced += fasta_seq[exon_start - 1 : exon_end - 1]
                else:
                    effective_length = (exon_end - exon_start + 1) + slippage
                    spliced += fasta_seq[exon_start - 1 : exon_start - 1 + effective_length]

        # Reverse complement to get mRNA orientation
        spliced = reverse_complement(spliced)
        return spliced


# =============================================================================
# CODON EXTRACTION - SIMPLIFIED (CDS is already in mRNA orientation)
# =============================================================================

def get_codon_and_position(cds_offset, spliced_cds, alt_base, strand):
    """
    Extract reference and alt codons from the spliced CDS.
    
    Since spliced_cds is already in mRNA orientation (5'->3' of transcript),
    we just extract the codon directly. No reverse complementing needed here.
    
    For reverse strand genes, the alt_base must already be complemented
    BEFORE calling this function.
    """
    codon_position = cds_offset % 3
    codon_start = cds_offset - codon_position

    ref_codon = spliced_cds[codon_start : codon_start + 3]

    alt_codon = list(ref_codon)
    alt_codon[codon_position] = alt_base
    alt_codon = "".join(alt_codon)

    protein_position = (codon_start // 3) + 1

    return ref_codon, alt_codon, codon_position, protein_position


# =============================================================================
# HELPER FUNCTIONS
# =============================================================================

def get_exon_number(variant_pos, cds_exons, strand):
    """
    Get the exon number in transcript order.
    For forward strand, exon 1 is the first genomic exon.
    For reverse strand, exon 1 is the LAST genomic exon.
    """
    if strand == "+":
        for idx, (exon_start, exon_end) in enumerate(cds_exons):
            if variant_pos >= exon_start and variant_pos <= exon_end:
                return idx + 1
    else:
        reversed_exons = list(reversed(cds_exons))
        for idx, (exon_start, exon_end) in enumerate(reversed_exons):
            if variant_pos >= exon_start and variant_pos <= exon_end:
                return idx + 1
    return 1


def is_valid_codon(codon):
    valid_bases = set("ACGTRYKMSWHBVDN")
    return len(codon) == 3 and all(c.upper() in valid_bases for c in codon)


# =============================================================================
# ARGUMENT PARSING AND FILE LOADING
# =============================================================================

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("avinput")
    parser.add_argument("GFF")
    parser.add_argument("FASTA")

    args = parser.parse_args()

    gff = pd.read_csv(args.GFF, sep="\t", header=None, skiprows=[i for i in range(0, 2)])

    variantFunction = pd.read_csv(args.avinput, sep="\t", header=None)

    variantFunctionName = args.avinput
    variantFunctionName = os.path.basename(variantFunctionName)

    variantFunctionName = re.sub(".fastq.avinput", "", str(variantFunctionName))

    # Load FASTA and convert to plain string for reliable slicing
    for fasta in SeqIO.parse(args.FASTA, "fasta"):
        print("")

    fasta_sequence = str(fasta.seq)

    variantFunction.insert(0, "", "")
    variantFunction.insert(0, " ", "")

    numRegions = 0

    for i in range(len(gff.index)):
        if pd.isnull(gff.iloc[i, 3]) != True:
            numRegions = numRegions + 1
        else:
            break

    test = 0
    df = None

    # =============================================================================
    # FIRST PASS: Annotate variants as "exonic" and assign gene names
    # =============================================================================

    for j in range(numRegions):
        if gff.iloc[j, 2] == "CDS":
            for k in range(len(variantFunction)):
                proteinName = gff.iloc[j, 8]
                proteinName = proteinName.replace("ID=CDS:", "gene:")
                sep = ";"
                proteinName = proteinName.split(sep, 1)[0]

                if variantFunction.iloc[k, 8] >= int(
                    gff.iloc[j, 3]
                ) and variantFunction.iloc[k, 8] <= int(gff.iloc[j, 4]):
                    if variantFunction.iloc[k, 0] == "":
                        variantFunction.iloc[k, 0] = "exonic"
                        variantFunction.iloc[k, 1] = proteinName

                    elif variantFunction.iloc[k, 0] != "":
                        if test == 0:
                            df = pd.DataFrame(variantFunction.loc[[k]])
                            df.iloc[0, 1] = proteinName

                        df2 = pd.DataFrame(variantFunction.loc[[k]])
                        df2.iloc[0, 1] = proteinName

                if df is not None:
                    df = df.iloc[1:, :]

    variantFunction = pd.concat([variantFunction, df])

    variantFunction.to_csv("varriantfunction.csv", index=False, header=False)
    variantFunction.to_csv(
        variantFunctionName + ".varriant_function",
        index=False,
        sep="\t",
        encoding="utf-8",
        header=False,
    )

    variantFunction.insert(0, "   ", "")
    variantFunction.insert(0, "    ", "")
    variantFunction.insert(19, "     ", "")

    # =============================================================================
    # SECOND PASS: Annotate mutation type, amino acid change, protein position
    # =============================================================================

    gene_cds_map, gene_strand_map, gene_coords_map = build_gene_cds_map(gff)

    for l in range(len(variantFunction)):
        if variantFunction.iloc[l, 2] != "exonic":
            continue

        variantFunction.iloc[l, 0] = "line" + str(l + 1)

        proteinName = variantFunction.iloc[l, 3]
        proteinName2 = proteinName.replace("gene:", "transcript:")

        gene_name = proteinName.replace("gene:", "").strip()

        if gene_name not in gene_strand_map:
            continue

        strand = gene_strand_map[gene_name]
        gene_start, gene_end = gene_coords_map.get(gene_name, (None, None))

        if gene_start is None:
            continue

        reverseComplement = "Yes" if strand == "-" else "No"
        variantFunction.iloc[l, 19] = reverseComplement

        cds_exons = gene_cds_map.get(gene_name, None)
        if cds_exons is None:
            cds_exons = [(gene_start, gene_end)]

        # --- Handle indels ---
        ref_allele = str(variantFunction.iloc[l, 12]).strip()
        alt_allele = str(variantFunction.iloc[l, 13]).strip()

        if len(ref_allele) > len(alt_allele):
            deleted_len = len(ref_allele) - len(alt_allele)
            if deleted_len % 3 == 0:
                variantFunction.iloc[l, 1] = "nonframeshift deletion"
            else:
                variantFunction.iloc[l, 1] = "frameshift deletion"

            # Add transcript annotation for indels
            variant_pos = int(variantFunction.iloc[l, 10])
            exon_number = get_exon_number(variant_pos, cds_exons, strand)
            variantFunction.iloc[l, 3] = (
                proteinName
                + ":"
                + proteinName2
                + ":exon"
                + str(exon_number)
                + ":c."
            )
            continue

        if len(ref_allele) < len(alt_allele):
            inserted_len = len(alt_allele) - len(ref_allele)
            if inserted_len % 3 == 0:
                variantFunction.iloc[l, 1] = "nonframeshift insertion"
            else:
                variantFunction.iloc[l, 1] = "frameshift insertion"

            # Add transcript annotation for indels
            variant_pos = int(variantFunction.iloc[l, 10])
            exon_number = get_exon_number(variant_pos, cds_exons, strand)
            variantFunction.iloc[l, 3] = (
                proteinName
                + ":"
                + proteinName2
                + ":exon"
                + str(exon_number)
                + ":c."
            )
            continue

        # --- Handle SNVs ---
        variant_pos = int(variantFunction.iloc[l, 10])

        ref_base_check = str(variantFunction.iloc[l, 12]).strip().upper()
        if ref_base_check == "N" or ref_base_check == "0":
            continue

        exon_number = get_exon_number(variant_pos, cds_exons, strand)

        cds_offset = get_cds_offset_with_slippage(variant_pos, cds_exons, strand)

        if cds_offset < 0:
            continue

        # Build spliced CDS in mRNA orientation
        spliced_cds = build_spliced_cds_with_slippage(fasta_sequence, cds_exons, strand)

        codon_position = cds_offset % 3
        codon_start = cds_offset - codon_position

        if codon_start + 3 > len(spliced_cds):
            continue

        # Handle "0" placeholders
        if "0" in str(variantFunction.iloc[l, 8]):
            variantFunction.iloc[l, 8] = str(variantFunction.iloc[l, 8]).replace(
                "0", str(variantFunction.iloc[l, 13])
            )
        if "0" in str(variantFunction.iloc[l, 7]):
            variantFunction.iloc[l, 7] = str(variantFunction.iloc[l, 7]).replace(
                "0", str(variantFunction.iloc[l, 12])
            )

        alt_base = str(variantFunction.iloc[l, 8]).strip().upper()

        if alt_base == "-":
            continue

        # KEY FIX: For reverse strand genes, complement the alt base
        # because VCF reports in genomic orientation but our spliced CDS
        # is in mRNA orientation (reverse complement of genomic)
        if strand == "-":
            alt_base = complement_base(alt_base)

        ref_codon_seq = spliced_cds[codon_start : codon_start + 3]
        if "-" in ref_codon_seq:
            continue

        # No strand argument needed - CDS is already in mRNA orientation
        ref_codon, alt_codon, codon_pos, protein_position = get_codon_and_position(
            cds_offset, spliced_cds, alt_base, strand
        )

        if not is_valid_codon(ref_codon):
            continue
        if not is_valid_codon(alt_codon):
            continue

        aminoNum = str(cds_offset + 1)
        proteinNum = str(protein_position)

        before = degenerative(ref_codon)
        after = degenerative(alt_codon)

        if len(before) == 3 and all(c in "ACGT" for c in before):
            before = translate(before)
        if len(after) == 3 and all(c in "ACGT" for c in after):
            after = translate(after)

        # For HGVS c. notation on reverse strand, complement the displayed bases
        hgvs_ref_base = str(variantFunction.iloc[l, 7]).strip()
        hgvs_alt_base = str(variantFunction.iloc[l, 8]).strip()
        if strand == "-":
            hgvs_ref_base = complement_base(hgvs_ref_base.upper())
            hgvs_alt_base = complement_base(hgvs_alt_base.upper())

        variantFunction.iloc[l, 3] = (
            proteinName
            + ":"
            + proteinName2
            + ":exon"
            + str(exon_number)
            + ":c."
            + hgvs_ref_base
            + aminoNum
            + hgvs_alt_base
            + ":p."
            + before
            + proteinNum
            + after
            + ","
        )

        if before == after:
            variantFunction.iloc[l, 1] = "synonymous SNV"
        elif after == "*" and before != "*":
            variantFunction.iloc[l, 1] = "stopgain"
        elif before == "*" and after != "*":
            variantFunction.iloc[l, 1] = "stoploss"
        else:
            variantFunction.iloc[l, 1] = "nonsynonymous SNV"


    # =============================================================================
    # OUTPUT
    # =============================================================================

    variantFunction = variantFunction[variantFunction.iloc[:, 2] == "exonic"]
    variantFunction.drop(variantFunction.columns[2], axis=1, inplace=True)

    variantFunctionName = re.sub(".avinput", "", str(variantFunctionName))
    variantFunction.to_csv("varriantfunction_exonic.csv", index=False, header=False)

    if ".fastq" in variantFunctionName:
        variantFunction.to_csv(
            variantFunctionName + ".exonic_variant_function",
            index=False,
            sep="\t",
            encoding="utf-8",
            header=False,
        )
    else:
        variantFunction.to_csv(
            variantFunctionName + ".fastq.exonic_variant_function",
            index=False,
            sep="\t",
            encoding="utf-8",
            header=False,
        )
