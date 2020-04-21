import glob
import os
import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo

def align_by_phylotype(input_fasta,input_cluster,input_metadata,output_folder):
    metadata_dic = {}
    phylotype_dic = {}
    seq_dic = {}
    consensus_dic = {}

    with open(input_metadata,"r") as f:
        reader = csv.DictReader(f)
        metadata = [r for r in reader]

    for items in metadata:
        metadata_dic[items["header"]] = items["lineage"]

    for record in SeqIO.parse(input_fasta, 'fasta'):
        seq_dic[record.id]= record.seq

    with open(input_cluster,"r") as f:
        for line in f:
            phylotype_dic[line.rstrip()] = []

    for cluster in phylotype_dic.keys():
        for seq_id,phylotype in metadata_dic.items():
            cluster_type = cluster.split(".")
            cluster_length = len(cluster_type)
            phylo_type = phylotype.split(".")
            if len(phylo_type) < cluster_length:
                continue
            if phylo_type[:cluster_length] == cluster_type:
                if seq_id in seq_dic.keys():
                    phylotype_dic[cluster].append([seq_id,seq_dic[seq_id],phylotype])
                    del seq_dic[seq_id]

    print("Clade","Number of Sequences")
    for key,value in phylotype_dic.items():
        print(key,len(value))

    log_file = open(output_folder+"align_phylo.log","w")
    for seq in seq_dic.keys():
        log_file.write("Sequence " + seq + " with lineage " + metadata_dic[seq] +
        " did not fall into any of the phylotypes stated in the cluster file.\n")

    for key in phylotype_dic.keys():
        if len(phylotype_dic[key]) > 2:
            outfile_name = output_folder + "lineage_" + key + ".fasta"
            outfile = open(outfile_name,"w")
            for sequences in phylotype_dic[key]:
                record = SeqRecord(sequences[1],id=sequences[0],description="")
                SeqIO.write(record, outfile, "fasta")
            outfile.close()
            alignment_name = outfile_name[:-6] + "_alignment.fasta"
            align_command = "mafft " + outfile_name + " > " + alignment_name
            os.system(align_command)
            os.remove(outfile_name)
            alignment = AlignIO.read(alignment_name, 'fasta')
            consensus_name = key + "_consensus"
            summary_align = AlignInfo.SummaryInfo(alignment)
            consensus_seq = summary_align.dumb_consensus(threshold=0.0,ambiguous='N')
            consensus_dic[consensus_name] = consensus_seq
        else:
            log_file.write("Phylotype " + key + "does not have 2 or more sequences for an alignment to work.")
    log_file.close()

    consensus_file = open(output_folder+"lineage_consensus.fasta","w")
    for key, value in consensus_dic.items():
        record = SeqRecord(value,id=key,description="")
        SeqIO.write(record, consensus_file, "fasta")
    consensus_file.close()
