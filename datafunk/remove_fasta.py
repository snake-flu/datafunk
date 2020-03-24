from Bio import SeqIO

def filter_list(input_filter):
    filter_dictionary = {}
    with open(input_filter,"r") as filter_file:
        for line in filter_file:
            comment_pos = line.find("#")
            if comment_pos == -1:
                filter_dictionary[line.rstrip()] = "No Comment"
            else:
                filter_dictionary[line[:comment_pos].rstrip()] = line[comment_pos+1:].rstrip()
    return filter_dictionary

def remove_fasta(input_fasta, filter_dictionary, output_file="filtered_file.fasta"):
    filtered_file = open(output_file,"w")
    log_file = open(output_file+".log","w")
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id not in filter_dictionary.keys():
            SeqIO.write(record, filtered_file, "fasta")
        else:
            log_file.write("Filtered Sequence: " + record.id + " due to " + filter_dictionary[record.id] + "\n")
    filtered_file.close()
    log_file.close()
                