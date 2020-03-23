from Bio import SeqIO

def make_record_dict(fasta):
    record_dict = {}
    for record in SeqIO.parse(fasta,"fasta"):
        iqtree_id = str(record.id).replace("|","_").replace("/","_") 
        record_dict[iqtree_id]=record.id
    return record_dict

def fix_names(fasta, tree, out):
    print("\n*** Fixing the taxa names ***\n")

    record_dict = make_record_dict(str(args.fasta))

    fw = open(str(out),"w")
    
    with open(str(tree),"r") as f:
        for l in f:
            l =l.rstrip("\n")
            new_l = l
            for record in record_dict:
                if record in new_l:
                    
                    new_l = new_l.replace(record, record_dict[record])
                else:
                    print("Missed this one:",record)
                    print("\n")

            fw.write(new_l + '\n')
    fw.close()




