from Bio import SeqIO
import sys, os


WH04_align = SeqIO.read(os.path.dirname(os.path.realpath(__file__)) + '/resources/WH04_aligned.fa', 'fasta')


def get_pairwise_difference(seq1, seq2):
    if len(seq1) != len(seq2):
        sys.exit('unequal sequence lengths')
    total=0
    diff=0
    for i in range(len(seq1)):
        a = seq1[i]
        b = seq2[i]
        if a in ['A','C', 'G', 'T'] and b in ['A','C', 'G', 'T']:
            total +=1
            if a != b:
                diff+=1
    return(total, diff, len(seq1))


def distance_per_genome(distance_info):
    """
    distance_info is a tuple : (comparisons, differences, genome_length)
    """
    comparisons = distance_info[0]
    differences = distance_info[1]
    genome_length = distance_info[2]

    if differences == 0:
        return(0)

    distance_per_genome = differences / comparisons * genome_length

    return(distance_per_genome)


def read_metadata(file, sep = ','):
    metadata = {}
    with open(file, 'r') as f:
        First = True
        for line in f:
            l = line.rstrip().split(sep)
            if First:
                header = l
                if not all(x in l for x in ['edin_omitted', 'subsample_omit', 'sequence_name', 'edin_epi_week']):
                    sys.exit('required columns not found in metadata')
                First = False
                continue
            d = {x:y for x,y in zip(header, l)}
            if d['edin_omitted'] == 'True':
                continue
            if d['subsample_omit'] == 'True':
                continue
            seq_name = d['sequence_name']
            if seq_name in metadata:
                warnings.warn('duplicate entry in ' + file + ', ignoring ' + seq_name)
                del metadata[seq_name]
                continue
            metadata[seq_name] = d

    return(metadata)


def distance_to_root(fasta_file, metadata_file):
    metadata = read_metadata(metadata_file)
    fasta = SeqIO.parse(fasta_file, 'fasta')

    out = open('distances.txt', 'w')
    out.write('sequence_name\tepi_week\tdistance_WH04\n')
    for record in fasta:
        id = record.id
        if id not in metadata:
            warnings.warn(id + ' not found in ' + metadata_file)

        metadata_line = metadata[id]

        seq = record.seq

        distance_info = get_pairwise_difference(WH04_align.seq, seq)
        distance = distance_per_genome(distance_info)

        try:
            epi_week = int(metadata_line['edin_epi_week'])
        except:
            epi_week = None

        if epi_week:
            out.write(id + '\t' + str(epi_week) + '\t' + str(round(distance,6)) + '\n')

    pass


















    #
