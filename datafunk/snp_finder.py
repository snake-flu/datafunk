import argparse
import collections
from Bio import AlignIO
import os
import csv

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
    reference = ''
    for record in aln:
        if "WH04" in record.id:
            reference = record
    return reference


def get_all_snps(alignment_file, snp, position, outfile, label_dict):
    aln = AlignIO.read(alignment_file, "fasta")

    reference = get_reference(aln)
    if type(reference) == str:
        sys.stderr.write('Error: couldnt find ref in file')
        sys.exit(-1)

    for record in aln:
        if record.id != reference.id:
            nucleotide = find_snp(reference.seq, record.seq, position)
            if nucleotide in label_dict:
                outfile.write(f"{record.id},{snp},{nucleotide},{label_dict[nucleotide]}\n")
            else:
                outfile.write(f"{record.id},{snp},{nucleotide},\n")


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
        fw.write("name,location,nucleotide,label\n")
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
                for k in label_dict:
                    print(f"Nuc: {k}, Label:{label_dict[k]}")
                location = int(row["location"])
                get_all_snps(alignment_file, row["name"], location, fw, label_dict)
