from collections import defaultdict
from collections import Counter
import os
import sys
import glob
import operator

class taxon():

    def __init__(self,name, lineage, acctrans):
        self.id = name
        self.lineage = lineage
        self.acctrans = acctrans


class lineage():

    def __init__(self, name, taxa):
        self.id = name
        self.taxa = taxa
        self.acctrans_designations = set()

        self.split = False
        self.merge = False


def make_taxon_objects(input_dir):

    taxon_list = []
    introduction_int_list = []
    intros_to_taxa = defaultdict(list)

    intro_acctrans = defaultdict(set)
    acctrans_to_intro = defaultdict(set)

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
                    intro_name = toks[3]
                    acctrans = toks[4]
                    new_taxon = taxon(seq_name, intro_name, acctrans)
                    taxon_list.append(new_taxon)

                    if intro_name != "":
                        intro_acctrans[intro_name].add(acctrans)
                        acctrans_to_intro[acctrans].add(intro_name)
                        intros_to_taxa[intro_name].append(new_taxon)
                        introduction_int_list.append(int(intro_name.lstrip("UK")))

    introduction_int_list = sorted(introduction_int_list)

    return introduction_int_list, taxon_list, intros_to_taxa, acctrans_to_intro, intro_acctrans


def deal_with_new_lineages(introduction_int_list, taxon_list, intros_to_taxa, acctrans_to_intro, intro_acctrans):
    
    new_acctrans_to_lineage = {}
    unclear_taxa = []

    for tax in taxon_list: #Assigning taxa to lineages that hasn't been done in the tree
        if tax.lineage== "":
                
                possible_intro = acctrans_to_intro[tax.acctrans]
                
                if len(possible_intro) == 0: #this will be if the introduction is genuinely new
                    if tax.acctrans not in new_acctrans_to_lineage.keys():

                        introduction_prep = (introduction_int_list[-1] + 1)
                        introduction = "UK" + str(introduction_prep)
                        introduction_int_list.append(introduction_prep)

                        new_acctrans_to_lineage[tax.acctrans] = introduction

                        tax.lineage = introduction
                        intros_to_taxa[introduction].append(tax)
                        intro_acctrans[introduction].add(tax.acctrans)
                    
                    else: #if it's not the sequence to be added to a new lineage
    
                        introduction = new_acctrans_to_lineage[tax.acctrans]

                        tax.lineage = introduction
                        
                        intros_to_taxa[introduction].append(tax)
                        intro_acctrans[introduction].add(tax.acctrans)

                elif len(possible_intro) > 1: #usually a merged or split UK lineage 
                    tax.unclear = True
                    unclear_taxa.append(tax)
                
                else: #if it just hasn't been labelled but it's clear it belongs to an existing lineage
                    introduction = list(possible_intro)[0]
                    
                    tax.lineage = introduction
                    intros_to_taxa[introduction].append(tax)
                    intro_acctrans[introduction].add(tax.acctrans)

    return intros_to_taxa, intro_acctrans, unclear_taxa

def make_lineage_objects(intros_to_taxa, intro_acctrans):

    lineage_objects = []
    lineage_dict = {}

    for intro, taxa in intros_to_taxa.items():
        
        i_o = lineage(intro, taxa)
        i_o.acctrans_designations = intro_acctrans[i_o.id]
        
        lineage_objects.append(i_o)

        # if len(i_o.acctrans_designations) > 1:
        #     i_o.split = True

        lineage_dict[intro] = i_o

    return lineage_objects, lineage_dict

def find_merged(lineage_objects):

    acctran_dict = defaultdict(list)
    acctran_to_merge = defaultdict(list)
    merged = {}
    merge_count = 0
    skip_list = []

    for lin in lineage_objects:         
        for i in lin.acctrans_designations: 
            acctran_dict[i].append(lin)
           
    for key, value in acctran_dict.items():
        if len(value) > 1:
            for i in value:
                merged[i.id] = True
                i.merge = True
                merge_count += 1
                if i not in acctran_to_merge[key]:
                    acctran_to_merge[key].append(i)

    return merged, acctran_dict, acctran_to_merge, merge_count

def deal_with_merged(acctran_to_merge, lineage_object_dict, unclear_taxa, lineage_objects, introduction_int_list):
    #selecting new names and making new lineage object and also updating taxon defs
    remove_list = set()
    already_used = set()
    taxa_already_assigned = []

    for acctran, lineages in acctran_to_merge.items():
        
        taxa = []
        numbers = []
        lineages_present = []
       
    
        for lin in lineages:
            for tax in lin.taxa:
                if tax.acctrans == acctran:
                    taxa.append(tax)
            
            remove_list.add(lin)
            numbers.append(int(lin.id.lstrip("UK")))
            
            
        for tax in taxa:
            lineages_present.append(tax.lineage)
        
        lineage_count = Counter(lineages_present)
        freq_order = lineage_count.most_common()
        
        #if two equal sized lineages
        possibles = []
        highest_count = freq_order[0][1]
        for i in freq_order:
            if i[1] == highest_count:
                number = int(i[0].lstrip("UK"))
                possibles.append(number)
        
        if len(possibles) > 1:
            new_name = "UK" + str(min(numbers))      
        else:
            new_name = freq_order[0][0]
          
        #Check it hasn't already been used
        count = 0
        new_name = lineage_count.most_common()[count][0]
        while new_name in already_used and count <= len(lineage_count)-1:
            count += 1
            new_name = lineage_count.most_common()[count][0]
        
        #If all of the lineage names have been used before
        if count > len(lineage_count)-1:
            introduction_prep = (introduction_int_list[-1] + 1)
            new_name = "UK" + str(introduction_prep)
            introduction_int_list.append(introduction_prep)
        
        
        already_used.add(new_name)

        test_taxa = set(taxa)
        if len(test_taxa) != len(taxa):
            print("error in merging" + str(lineages) + " " + str(len(test_taxa)) + " " + str(len(taxa)))


        for tax in unclear_taxa:
            if tax.acctrans == acctran:
                tax.lineage= new_name
                taxa.append(tax)
                
        taxa_already_assigned.extend(taxa)

        i_o = lineage(new_name, taxa)
        i_o.acctrans_designations.add(acctran)

        for tax in i_o.taxa:
            tax.lineage = i_o.id

        lineage_object_dict[i_o.id] = i_o      

        lineage_objects.append(i_o)
       
    return lineage_object_dict, lineage_objects, introduction_int_list, taxa_already_assigned, remove_list

def deal_with_split(lineage_objects, lineage_dict, acctran_dict, taxa_already_assigned, introduction_int_list, remove_list):
    
    acctran_for_split = defaultdict(list)
    lin_to_acctran = defaultdict(list)
    lin_to_size = defaultdict(dict)
    acc_to_count = {}
    max_acc_to_lin = {}
    
    lin_to_biggest_acctran = {}
    
    relevant_lineages = set()
    
    for lin in lineage_objects:
        if len(lin.acctrans_designations) > 1:
            lin.split = True
         
    
        
    for acc, lins in acctran_dict.items():
        if len(lins) == 1 and lins[0].split:
            remove_list.add(lins[0])
            relevant_lineages.add(lins[0])
            for tax in lins[0].taxa:
                if tax not in taxa_already_assigned:
                    acctran_for_split[tax.acctrans].append(tax)
                    taxa_already_assigned.append(tax)
                    
 
            acc_to_count[acc] = len(acctran_for_split[acc])
             
    for lin in relevant_lineages:
        acc_dict = {}
        for acc in lin.acctrans_designations:
            if acc in acc_to_count.keys():
                acc_dict[acc] = acc_to_count[acc]
        lin_to_size[lin] = acc_dict
        
        max_acc = max(acc_dict.items(), key=operator.itemgetter(1))[0]
        
        lin_to_biggest_acctran[lin] = max_acc
        max_acc_to_lin[max_acc] = lin.id
                    
    for acc, taxa in acctran_for_split.items():
        
        if acc in max_acc_to_lin.keys():
            new_name = max_acc_to_lin[acc]
        else:
            introduction_prep = (introduction_int_list[-1] + 1)
            new_name = "UK" + str(introduction_prep)
            introduction_int_list.append(introduction_prep)
        
        i_o = lineage(new_name, taxa)
        i_o.acctrans_designations.add(acc)
        
        lineage_objects.append(i_o)
        lineage_dict[new_name] = i_o
        
        for tax in i_o.taxa:
            tax.lineage = i_o.id
    
    return lineage_objects, lineage_dict, remove_list, acctran_for_split


def write_to_file(lineage_objects, outfile):

    fw = open(outfile, 'w')
    fw.write("taxon,uk_lineage\n")
    for lin in lineage_objects:
        for tax in lin.taxa:
            fw.write(tax.id + "," + tax.lineage + "\n")

    fw.close()


def curate_lineages(input_dir, outfile):

    introduction_int_list, taxon_list, intros_to_taxa, acctrans_to_intro, intro_acctrans = make_taxon_objects(input_dir)

    intros_to_taxa, intro_acctrans, unclear_taxa = deal_with_new_lineages(introduction_int_list, taxon_list, intros_to_taxa, acctrans_to_intro, intro_acctrans)

    lineage_objects, lineage_dict = make_lineage_objects(intros_to_taxa, intro_acctrans)

    merged, acctran_dict, acctran_to_merge, merge_count = find_merged(lineage_objects)

    lineage_object_dict, lineage_objects, introduction_int_list, taxa_already_assigned, remove_list = deal_with_merged(acctran_to_merge, lineage_dict, unclear_taxa, lineage_objects, introduction_int_list)

    lineage_objects, lineage_dict, remove_list, acctran_for_split = deal_with_split(lineage_objects, lineage_dict, acctran_dict, taxa_already_assigned, introduction_int_list, remove_list)

    for i in remove_list:
        if i in lineage_objects:
            lineage_objects.remove(i)

    write_to_file(lineage_objects, outfile)
