from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

# import os
import sys


# def parse_gtf(file, attribute_key_value_separator = '='):
#     sep = attribute_key_value_separator
#     d = []
#     with open(file, 'r') as f:
#         for line in f:
#             if line[0] == '#':
#                 continue
#             l = line.rstrip()
#             if len(l) == 0:
#                 continue
#             l = l.split('\t')
#
#             fields = ['seqname', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attribute']
#             values = l[0:8] + [{x[0]:x[1] for x in [(y.split(sep)[0],y.split(sep)[1]) for y in l[8].split(';')]}]
#
#             d.append({x: y for x,y in zip(fields, values)})
#
#     return(d)


# def get_one_feature(seq, start, end):
#     """
#     start and end are 1-based inclusive cordinates from a gtf file
#     """
#     return(seq[start-1:end])


# def get_all_CDSs(sequence, gtf_list):
#     AA_list = []
#     for entry in gtf_list:
#         if entry['feature'] == 'CDS':
#             start = int(entry['start'])
#             end = int(entry['end'])
#             s = get_one_feature(seq = sequence, start = start, end = end)
#             AA = get_AA(seq = s)
#             AA_list.append(AA)
#
#     return(AA_list)





CDS_coordinates = [((266, 13468), (13468, 21555)),
                    (21563,	25384),
                    (25393,	26220),
                    (26245,	26472),
                    (26523,	27191),
                    (27202,	27387),
                    (27394,	27759),
                    (27894,	28259),
                    (28274,	29533),
                    (29558,	29674)]


def get_CDS(fasta_in, fasta_out, translate = False):
    # gtf = os.path.dirname(os.path.realpath(__file__)) + '/resources/mn908947.3.gff3'
    # gtf_info = parse_gtf(gtf)


    if fasta_out:
        out = open(fasta_out, 'w')
    else:
        out = sys.stdout

    # # test stuff with reference
    # ref_aa_file = SeqIO.parse('/Users/ben/Dropbox/Biology_Edinburgh_2/Rambaut/cds/mn908947.3.aa', 'fasta')
    # ref_aa = []
    # for record in ref_aa_file:
    #     ref_aa.append(str(record.seq))
    # out.write('>reference\n')
    # out.write(''.join(ref_aa) + '\n')


    fasta = SeqIO.parse(fasta_in, 'fasta')

    for record in fasta:
        seqs = []
        for coords in CDS_coordinates:
            if isinstance(coords[0], tuple):
                seq = record.seq[coords[0][0] - 1:coords[0][1]] + record.seq[coords[1][0] - 1:coords[1][1]]
            elif isinstance(coords[0], int):
                seq = record.seq[coords[0] - 1:coords[1]]
            seqs.append(seq)

        nostops = [seq[:-3] for seq in seqs]

        if translate:
            out.write(record.id + '\n')
            out.write(''.join([str(seq.translate()) for seq in nostops]) + '\n')

        else:
            out.write(record.id + '\n')
            out.write(''.join([str(seq) for seq in nostops]) + '\n')

    if fasta_out:
        out.close()

    pass
