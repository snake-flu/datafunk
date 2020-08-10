from Bio import SeqIO
import sys, random

def bootstrap(fasta_in, n = 1, output_prefix = "boostrap_"):
    temp_alignment = SeqIO.parse(fasta_in, "fasta")
    temp_record = next(temp_alignment)

    # l is the number of sites in the alignment
    l = len(temp_record.seq)

    # each i is a bootstrap
    for i in range(n):
        alignment = SeqIO.parse(fasta_in, "fasta")

        # for each bootstrap, make one sample vector, the same size as the alignment is wide,
        # by drawing sites with replacement
        V = random.choices(range(l), k=l)

        # then iterate over each entry in the fasta file and rewrite a
        # bootstrapped sequence (using the same indices each time)
        with open(output_prefix + str(i + 1) + ".fasta", "w") as f_out:
            for record in alignment:
                id = record.id
                seq = record.seq
                new_seq = ''.join([seq[x] for x in V])
                f_out.write(">" + id + "\n")
                f_out.write(new_seq + "\n")
    pass
