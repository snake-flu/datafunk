import csv
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def import_file(input_metadata,input_fasta):
    UK_metadata = {}
    UK_all = {}
    UK_seq = {}
    with open(input_metadata,"r") as f:
        reader = csv.DictReader(f)
        metadata = [r for r in reader]
        
    for items in metadata:
        if items["UK"] == "True":
            seq_id = items["taxon"][:items["taxon"].find("|")]
            UK_metadata[seq_id] = [items["taxon"],items["phylotype"],items["ACCTRAN"]]

    for record in SeqIO.parse(input_fasta, 'fasta'):
        seq_id = record.id[:record.id.find("|")].replace("/","_")
        if seq_id in list(UK_metadata.keys()):
            UK_seq[UK_metadata[seq_id][0]] = record.seq
        
    for value in UK_metadata.values():
        seq_id = value[0]
        phylotype = value[1]
        intro = value[2]
        if intro not in list(UK_all.keys()):
            UK_all[intro] = []
            UK_all[intro].append([seq_id,UK_seq[seq_id],phylotype])
        else:
            UK_all[intro].append([seq_id,UK_seq[seq_id],phylotype])   
    return UK_all

def split_phylotype(output_folder,phylotype_dictionary,clade_threshold=2):
    log_file = open(output_folder+"phylotype.log","w")
    for keys in phylotype_dictionary.keys():
        if len(phylotype_dictionary[keys]) >= clade_threshold:
            outfile = open(output_folder+str(keys)+".fasta","w")
            for seq in phylotype_dictionary[keys]:
                record = SeqRecord(seq[1],id=seq[0],description="")
                SeqIO.write(record, outfile, "fasta")
                log_file.write("Sequence "+seq[0]+" has been writen to phylotype "+keys+ " file \n")
        else:
            log_file.write("Phylotype "+keys+" Contain atleast "+str(clade_threshold)+ " sequences and not written to file for downstream analysis\n")
    log_file.close()
