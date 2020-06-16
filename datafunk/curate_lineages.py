from collections import defaultdict
from collections import Counter
import os
import sys
import glob
import operator

class taxon():

    def __init__(self,name, lineage, acctrans):
        self.id = name
        self.old_lineage = lineage
        self.acctrans = acctrans


def make_taxon_objects(input_dir):

    taxon_list = []
    acc_list = set()

    if input_dir.endswith("traits.csv"):
        list_traits_files = [input_dir]
    else:
        in_dir = os.path.dirname(input_dir)
        list_traits_files = glob.glob('%s/**/*traits.csv' %in_dir, recursive=True) \
                            + glob.glob('%s/*traits.csv' %in_dir, recursive=True)
    if len(list_traits_files) == 0:
        sys.exit("Found no traits files!!")
    print("Found traits files:", list_traits_files)

    for input_file in list_traits_files:
        with open(input_file) as f:
            next(f)
            for l in f:
                toks = l.strip("\n").split(",")
                if toks[1] == "UK": #this is in the traits.csv
                    seq_name = toks[0]
                    intro_name = toks[2]
                    acctrans = toks[6] #This is actually now deltrans, but easier to keep the names the same
                    new_taxon = taxon(seq_name, intro_name, acctrans)
                    taxon_list.append(new_taxon)

                    acc_list.add(acctrans)

    acc_to_tax = defaultdict(list)

    for tax in taxon_list:
        acc_to_tax[tax.acctrans].append(tax)

    acc_to_poss_lins = defaultdict(list)
    lins_to_acc = defaultdict(list)

    for k,v in acc_to_tax.items():
        for tax in v:
            acc_to_poss_lins[k].append(tax.old_lineage)
            
    for k,v in acc_to_poss_lins.items():
        for lineage in v:
            lins_to_acc[lineage].append(k)

    acc_name_counts = defaultdict(list)
    lin_acc_counts = defaultdict(list)

    for k,v in acc_to_poss_lins.items():
        new = Counter(v)
        acc_name_counts[k] = new

        
    for k,v in lins_to_acc.items():
        new = Counter(v)
        lin_acc_counts[k] = new  



    return taxon_list, acc_list, acc_name_counts, lin_acc_counts, acc_to_tax

def rename_lineages(acc_name_counts, lin_acc_counts):

    acc_final_name_dict = defaultdict(list) #this will all go into one, just for testing

    new_names = []

    for acc, poss_lins in acc_name_counts.items():
        
        relevant_lin = next(iter(poss_lins.items()))[0]
        
        if len(poss_lins) == 1 and relevant_lin != "" and len(lin_acc_counts[relevant_lin]) == 1:
                
            acc_final_name_dict[acc].append(relevant_lin)
        
        #For new lineages, so they have no old lineage to fall back on
        elif len(poss_lins) == 1 and relevant_lin == "":
            acc_final_name_dict[acc].append(relevant_lin) 
        else:
            new_name = None 

            for i in range(1, len(poss_lins)+1):
                more_deserving = False
                if new_name == None:
                    other_deserving = {}

                    query_lin = poss_lins.most_common(i)[i-1][0]
                    count = poss_lins.most_common(i)[i-1][1]

                    if query_lin == "": #so skip any empty strings here
                        continue

                    else:

                        for inner_acc, inner_poss_lins in acc_name_counts.items():

                            test_against_other = inner_poss_lins.most_common(1)[0][0]
                                        
                            if query_lin == test_against_other and acc != inner_acc:

                                other_deserving[inner_acc] = inner_poss_lins.most_common(1)[0][1]                          

                        if len(other_deserving) == 0:
                            new_name = query_lin

                        else:
                            for other_acc, other_count in other_deserving.items():
                                if other_count > count:
                                    more_deserving = True
                                    count = other_count


                        if more_deserving:
                            continue

                        else:
                            new_name = query_lin
                            
                    
            if new_name == None: #if we loop to the end of the thing and all the names are taken
                new_name = ""
                
            acc_final_name_dict[acc].append(new_name)
            
            
            new_names.append(new_name)

    return acc_final_name_dict, new_names

def deal_with_issues(acc_final_name_dict, new_names, lin_acc_counts, acc_name_counts):

    name_counter = Counter(new_names)
    problem_lins = []

    for name, count in name_counter.items():
        if count > 1 and name != "":
            problem_lins.append(name)

    problem_lin = defaultdict(list)

    for lin in problem_lins:
        for acc, count in acc_name_counts.items():
            if lin in count:
                problem_lin[lin] = lin_acc_counts[lin]
        
        
    for lin, acc_options in problem_lin.items():
        
        for acc in acc_options.keys():
            acc_final_name_dict[acc] = []

        winner = acc_options.most_common(1)[0][0]
        
        acc_final_name_dict[winner].append(lin)

        
        for other_acc in acc_options:
            if other_acc != winner:
                for other_options in acc_name_counts[other_acc]:
                    if other_options != "" and other_options not in new_names:
                        new_name = other_options
                    else:
                        new_name = ""
                    
                acc_final_name_dict[other_acc].append(new_name)

    return acc_final_name_dict

def name_new_lineages(acc_final_name_dict):

    used_names = []
    for key, value in acc_final_name_dict.items():
        if value[0] != "":
             used_names.append(int(value[0].lstrip("UK")))
                
    sorted_names = (sorted(used_names))
    test_counter = Counter(sorted_names)

    full_list = []
    usable_names = []

    for i in range(1,len(used_names)+1):
        full_list.append(i)
        
    for i in full_list:
        if i not in used_names:
            usable_names.append(i)

    acc_final = {}

    for acc, lin in acc_final_name_dict.items():
        if lin[0] == "":
            if len(usable_names) > 0:
                new_name_prep = usable_names[0]
                usable_names.remove(new_name_prep)
            else:
                new_name_prep = len(full_list) + 1
                full_list.append(new_name_prep)
          
            new_name = "UK" + str(new_name_prep)
            acc_final[acc] = new_name
        else:
            acc_final[acc] = lin[0]


    return acc_final, test_counter



def write_to_file(acc_to_tax, acc_final, outfile):

    new_tax_list = []
    for acc, tax_list in acc_to_tax.items():
        for tax in tax_list:
            tax.uk_lineage = acc_final[acc]
            new_tax_list.append(tax)

    for_counting_lin_size = defaultdict(list)

    for tax in new_tax_list:
        for_counting_lin_size[tax.uk_lineage].append(tax)

    lin_counts = {}
    for i,v in for_counting_lin_size.items():
        lin_counts[i] = len(v)

    counter_lin_counts = Counter(lin_counts)

    top_20_prep = counter_lin_counts.most_common(20)
    top_20 = []
    for i in top_20_prep:
        top_20.append(i[0])
    

    fw = open(outfile, 'w')
    fw.write("taxon,uk_lineage,acctrans,microreact_lineage\n")
    for tax in new_tax_list:
        if tax.uk_lineage in top_20:
            new_line = tax.id + "," + tax.uk_lineage + "," + tax.acctrans + "," + tax.uk_lineage + "\n"
        else:
            new_line = tax.id + ","  + tax.uk_lineage + "," + tax.acctrans + ",other\n"

        fw.write(new_line)

    fw.close()


def curate_lineages(input_dir, outfile):

    taxon_list, acc_list, acc_name_counts, lin_acc_counts, acc_to_tax = make_taxon_objects(input_dir)

    acc_final_name_dict, new_names = rename_lineages(acc_name_counts, lin_acc_counts)

    acc_final_name_dict = deal_with_issues(acc_final_name_dict, new_names, lin_acc_counts, acc_name_counts)

    acc_final, test_counter = name_new_lineages(acc_final_name_dict)

    for k,v in test_counter.items(): #Tests if a UK lineage name has been used more than once
        assert v == 1
    for acc, lin in acc_final_name_dict.items(): #Tests if a UK lineage has more than one acctrans assigned to it
        assert len(lin) == 1
        
    write_to_file(acc_to_tax, acc_final, outfile)

