from Bio import SeqIO
import os, sys

"""
The reference sequence, for genotyping:
"""
WuhanHu1 = SeqIO.read(os.path.dirname(os.path.realpath(__file__)) + '/resources/Wuhan-Hu-1.fa', 'fasta')


def parse_del_file(file):
    """
    input is in the format:
    start (1-based), length of deletion
    e.g.:
    1605,3

    l is a list of length-2 tuples with the format (position, length)

    it has the same number of entries as lines in file
    """

    ls = []

    with open(file, 'r') as f:
        for line in f:
            l = line.rstrip().split(',')
            pos, len = l

            ls = ls + [(int(pos), int(len))]

    return(ls)


def del_finder(fasta_in, fasta_out, del_file, genotypes_file, append_snp = False):
    """
    For every record in the query fasta file, for every deletion defined in del_file,
    genotype the record and write the genotype to a csv file. Optionally add the deletion/deletions'
    genotypes as a SNP at the end of the record in fasta_out
    """

    if fasta_out:
        f_out = open(fasta_out, 'w')

    g_out = open(genotypes_file, 'w')

    dels = parse_del_file(del_file)

    g_out.write("sequence_name," + ",".join(["del_" + str(x[0]) + "_" + str(x[1]) for x in dels]) + '\n')

    input = SeqIO.parse(fasta_in, 'fasta')

    for record in input:
        ID = record.id
        seq = str(record.seq).upper()

        if len(seq) != len(WuhanHu1):
            sys.exit("reference and query sequences are not the same length!")

        for entry in dels:
            pos = entry[0]
            len = entry[1]

            REF_allele = str(WuhanHu1.seq).upper()[pos - 1: pos - 1 + len]

            if seq[pos - 1: pos - 1 + len] == '-' * len:
                nuc = 'C'
                genotype = 'del'
            elif seq[pos - 1: pos - 1 + len] == REF_allele:
                nuc = 'A'
                genotype = 'ref'
            else:
                nuc = 'N'
                genotype = 'X'

            if append_snp:
                seq = seq + nuc

            # print(seq[pos - 1: pos - 1 + len], genotype, nuc)

        if fasta_out:
            f_out.write('>' + ID + '\n')
            f_out.write(seq + '\n')

        g_out.write(ID + "," + genotype + "\n")

    if fasta_out:
        f_out.close()

    g_out.close()
