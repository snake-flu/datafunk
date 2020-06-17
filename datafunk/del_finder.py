from Bio import SeqIO


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

    f_out = open(fasta_out, 'w')
    g_out = open(genotypes_file, 'w')

    dels = parse_del_file(del_file)

    g_out.write("sequence_name," + ",".join(["del_" + str(x[0]) + "_" + str(x[1]) for x in dels]) + '\n')

    input = SeqIO.parse(fasta_in, 'fasta')

    for record in input:
        ID = record.id
        seq = str(record.seq).upper()

        for entry in dels:
            pos = entry[0]
            len = entry[1]

            if seq[pos - 1: pos - 1 + len] == '-' * len:
                nuc = 'C'
                genotype = 'del'
            else:
                nuc = 'A'
                genotype = 'ref'

            if append_snp:
                seq = seq + nuc


            print(seq[pos - 1: pos - 1 + len], genotype, nuc)


        f_out.write('>' + ID + '\n')
        f_out.write(seq + '\n')

        g_out.write(ID + "," + genotype + "\n")

    f_out.close()
    g_out.close()
