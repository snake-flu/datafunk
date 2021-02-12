from Bio import SeqIO
import numpy as np
import sys, os


WH04_align = SeqIO.read(os.path.dirname(os.path.realpath(__file__)) + '/resources/WH04_aligned.fa', 'fasta')


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_pairwise_difference(seq1, seq2):
    if len(seq1) != len(seq2):
        sys.exit('unequal sequence lengths')
    total=0
    diff=0
    for i in range(len(seq1)):
        a = seq1[i].upper()
        b = seq2[i].upper()
        if a in ['A', 'C', 'G', 'T'] and b in ['A', 'C', 'G', 'T']:
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
                if not all(x in l for x in ['sequence_name', 'edin_epi_week']):
                    sys.exit('required columns not found in metadata')
                First = False
                continue
            d = {x:y for x,y in zip(header, l)}
            if 'edin_omitted' in header and d['edin_omitted'] == 'True':
                continue
            if 'subsample_omit' in header and d['subsample_omit'] == 'True':
                continue
            if 'why_excluded' in header and d['why_excluded']:
                continue
            seq_name = d['sequence_name']
            if seq_name in metadata:
                eprint('duplicate entry in ' + file + ', ignoring ' + seq_name)
                del metadata[seq_name]
                continue
            metadata[seq_name] = d

    return(metadata)


def get_epi_week_distance_stats(metadata):
    epi_weeks = {}
    for key in metadata:
        entry = metadata[key]
        try:
            epi_week = str(int(float(entry['edin_epi_week'])))
            distance = entry['distance']
        except:
            continue

        if epi_week not in epi_weeks:
            epi_weeks[epi_week] = {'distance': [distance]}
        else:
            epi_weeks[epi_week]['distance'] = epi_weeks[epi_week]['distance'] + [distance]

    for key in epi_weeks:
        mean = np.mean(epi_weeks[key]['distance'])
        std = np.std(epi_weeks[key]['distance'])

        epi_weeks[key]['mean'] = mean
        epi_weeks[key]['std'] = std

    return(epi_weeks)


def distance_to_root(fasta_file, metadata_file):
    metadata = read_metadata(metadata_file)
    fasta = SeqIO.parse(fasta_file, 'fasta')

    for record in fasta:
        id = record.id
        if id not in metadata:
            eprint(id + ' not found in ' + metadata_file)
            continue

        metadata_line = metadata[id]

        seq = record.seq
        print(len(seq), len(WH04_align.seq))

        distance_info = get_pairwise_difference(WH04_align.seq, seq)
        distance = distance_per_genome(distance_info)

        metadata[id]['distance'] = distance


    stats = get_epi_week_distance_stats(metadata)

    fasta = SeqIO.parse(fasta_file, 'fasta')
    out = open('distances.tsv', 'w')
    out.write('sequence_name\tepi_week\tepi_week_mean_distance\tepi_week_stdev_distance\tsample_distance\tdistance_stdevs\n')
    for record in fasta:
        id = record.id
        if id not in metadata:
            continue

        metadata_line = metadata[id]
        try:
            epi_week = str(int(float(metadata_line['edin_epi_week'])))
            dist = metadata[id]['distance']
        except:
            continue

        epi_week_mean_dist = stats[epi_week]['mean']
        epi_week_std_dist = stats[epi_week]['std']

        distance_std_units = (dist - epi_week_mean_dist) / epi_week_std_dist

        out.write(record.id + '\t' + epi_week + '\t' + str(round(epi_week_mean_dist, 4)) + '\t' + str(round(epi_week_std_dist, 4)) + '\t' + str(round(dist, 4)) + '\t' + str(round(distance_std_units, 4)) + '\n')

    out.close()

    pass


















    #
