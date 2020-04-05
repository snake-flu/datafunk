from Bio import SeqIO

def sequence_has_low_coverage(sequence, coverage_threshold):
    unaligned_seq = sequence.replace("-", "")
    high_covg_seq = unaligned_seq.replace("N", "")
    if len(unaligned_seq) == 0 or float(len(high_covg_seq)) / len(unaligned_seq) < coverage_threshold / 100.0:
        return True
    return False

def sequence_too_short(sequence, length_threshold):
    unaligned_seq = sequence.replace("-", "")
    if len(unaligned_seq) < length_threshold:
        return True
    return False

def filter_sequences(inpath, outpath, min_covg=None, min_length=None):
    print("\n*** Filtering out sequences with less than %i coverage ***\n" %min_covg)

    record_dict = SeqIO.index(inpath, "fasta")

    if outpath is None:
        outpath = inpath.replace(".fa",".filtered.fa")

    low_covg_seqs = []
    short_seqs = []

    with open(outpath, "w") as out_handle:
        for seq_name in record_dict:
            record_seq = str(record_dict[seq_name].seq)
            if min_covg and sequence_has_low_coverage(record_seq, min_covg):
                low_covg_seqs.append(seq_name)
                continue
            if min_length and sequence_too_short(record_seq, min_length):
                low_covg_seqs.append(seq_name)
                continue
            SeqIO.write(record_dict[seq_name], out_handle, "fasta")

        if min_covg:
            print("#Low coverage sequences:")
            for seq_name in low_covg_seqs:
                print(seq_name)
        if min_length:
            print("#Short/truncated sequences:")
            for seq_name in lshort_seqs:
                print(seq_name)
