from Bio import SeqIO
from Bio.Seq import Seq
# from Bio.Alphabet import generic_dna
import os, sys

"""
The reference sequence, for genotyping:
"""
WuhanHu1 = SeqIO.read(os.path.dirname(os.path.realpath(__file__)) + '/resources/Wuhan-Hu-1.fa', 'fasta')


def parse_AA_file(file):
    """
    input is in the format:
    start (1-based)
    e.g.:
    D614G,1605

    ls is a list of length-2 tuples with the format (name, position)
    position is the starting position of the codon in Wuhan-Hu-1 coordinates

    it has the same number of entries as lines in file
    """

    ls = []

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip().split(",")
            name, pos = l

            ls = ls + [(name, int(pos))]

    return(ls)


def AA_finder(fasta_in, AA_file, genotypes_file):
    """
    For every record in the query fasta file, for every codon start defined in AA_file,
    genotype the record's AA at that codon and write the genotype to a csv file.
    """

    g_out = open(genotypes_file, 'w')

    AAs = parse_AA_file(AA_file)

    g_out.write("sequence_name," + ",".join([x[0] for x in AAs]) + '\n')

    input = SeqIO.parse(fasta_in, 'fasta')

    for record in input:
        ID = record.id
        seq = record.seq

        if len(seq) != len(str(WuhanHu1.seq)):
            sys.exit("reference and query sequences are not the same length!")

        # a list of genotypes in case there is more than one allele in the file
        genotypes = []

        for entry in AAs:
            pos = entry[1]

            QUERY_seq = seq[pos - 1: pos + 2]

            if any([x in ['-', '?'] for x in QUERY_seq]):
                QUERY_allele = 'X'
            else:
                QUERY_allele = str(seq[pos - 1: pos + 2].translate())

            genotypes.append(QUERY_allele)

        g_out.write(ID + "," + ",".join(genotypes) + "\n")

    g_out.close()
