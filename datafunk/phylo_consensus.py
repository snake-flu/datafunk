from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio import AlignIO
from Bio.Align import AlignInfo
import glob

def create_consensus(input_folder,output_file="consensus_file.fasta"):
    UK_consensus={}
    for files in glob.glob(input_folder+"*.fasta"):
        alignment = AlignIO.read(files, 'fasta')
        slash_pos = files.rfind("/")
        consensus_name = files[slash_pos+1:-6]+"_consensus"
        summary_align = AlignInfo.SummaryInfo(alignment)
        consensus_seq = summary_align.gap_consensus(threshold=0.0,ambiguous='N')
        UK_consensus[consensus_name] = consensus_seq

    f = open(output_file,"w")
    for key, value in UK_consensus.items():
        record = SeqRecord(value,id=key,description="")
        SeqIO.write(record, f, "fasta") 
    f.close()