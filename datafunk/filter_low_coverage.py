from Bio import SeqIO

def remove_low_coverage_sequences(inpath, outpath, percent_threshold):
    print("\n*** Filtering out sequences with less than %i coverage ***\n" %percent_threshold)

    record_dict = SeqIO.index(inpath, "fasta")

    if outpath is None:
        outpath = inpath.replace(".fa",".filtered.fa")

    with open(outpath, "w") as out_handle:
        for seq_name in record_dict:
            record_seq = str(record_dict[seq_name].seq)
            unaligned_seq = record_seq.replace("-","")
            high_covg_seq = unaligned_seq.replace("N","")
            if len(unaligned_seq)==0 or float(len(high_covg_seq))/len(unaligned_seq) < percent_threshold/100.0:
                print("Sequence %s had (low) coverage %f" % (seq_name, float(len(high_covg_seq)) / len(unaligned_seq)))
                continue
            SeqIO.write(record_dict[seq_name], out_handle, "fasta")