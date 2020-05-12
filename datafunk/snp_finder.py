import argparse
import collections
from Bio import AlignIO
from Bio import SeqIO
import os
import csv
import sys

cwd = os.getcwd()

def find_snp(ref, member, position):
    index = 0
    snp = ""
    for i in range(len(ref)):
        if ref[i] != '-':
            index += 1

        if index == position:
            col = [ref[i], member[i]]
            snp = f"{col[1].upper()}"
    return snp


def get_reference(aln):
    reference = None
    for record in aln:
        if "WH04" in record.id:
            reference = record
    if reference is None:
        reference = SeqIO.read(os.path.dirname(os.path.realpath(__file__)) + '/resources/WH04_aligned.fa', 'fasta')
    if aln.get_alignment_length() != len(reference.seq):
        sys.stderr.write('Error: reference length is different from alignment length - have you trimmed already?!')
        sys.exit(-1)
    return reference


def get_all_snps(alignment_file, snp, position, outfile, label_dict):
    aln = AlignIO.read(alignment_file, "fasta")

    reference = get_reference(aln)
    if reference is None:
        sys.stderr.write('Error: couldnt find ref in file')
        sys.exit(-1)
    snp_dict = {}
    for record in aln:
        if record.id != reference.id:
            nucleotide = find_snp(reference.seq, record.seq, position)
            if nucleotide in label_dict:
                snp_dict[record.id] = label_dict[nucleotide]
            else:
                snp_dict[record.id] = "X"
    return snp_dict


def read_alignment_and_get_snps(alignment, snp_csv, outfile):

    alignment_file = os.path.join(cwd, alignment)
    if not os.path.exists(alignment_file):
        sys.stderr.write('Error: cannot find alignment file at {}\n'.format(alignment_file))
        sys.exit(-1)
    else:
        print(f"Reading in alignment file {alignment_file}.")

    snp_file = os.path.join(cwd, snp_csv)
    if not os.path.exists(snp_file):
        sys.stderr.write('Error: cannot find snp file at {}\n'.format(snp_file))
        sys.exit(-1)
    else:
        print(f"Reading in snp file {snp_file}.")

    with open(outfile, "w") as fw:

        tax_dict = collections.defaultdict(list)
        header = "name,"

        with open(snp_file, newline="") as csvfile:
            """
            name,location,nuc1,label1,nuc2,label2
            D614G,23403,A,G,D,G
            """
            reader = csv.DictReader(csvfile)
            for row in reader:

                label_dict = {
                    row["nuc1"]: row["label1"],
                    row["nuc2"]: row["label2"]
                }
                header += row["name"] + ','

                for k in label_dict:
                    print(f"Nuc: {k}, Label:{label_dict[k]}")

                location = int(row["location"])
                snp_dict = get_all_snps(alignment_file, row["name"], location, fw, label_dict)
                for record in snp_dict:
                    tax_dict[record].append(snp_dict[record])

            header = header.rstrip(',')
            fw.write(f"{header}\n")
            for record in tax_dict:
                snps = ",".join(tax_dict[record])
                line = f"{record},{snps}\n"
                fw.write(line)
