from collections import defaultdict
import os
import glob

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


def make_taxon_objects(input_dir):

    taxon_list = []
    introduction_int_list = []
    intros_to_taxa = defaultdict(list)

    intro_acctrans = defaultdict(set)
    acctrans_to_intro = defaultdict(set)

    list_traits_files = glob.glob('%s/**/*traits.csv' %input_dir, recursive=True)
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

    lineage_objects = []
    lineage_dict = {}

    for intro, taxa in intros_to_taxa.items():
        
        i_o = lineage(intro, taxa)
        i_o.acctrans_designations = intro_acctrans[i_o.id]
        
        lineage_objects.append(i_o)

        if len(i_o.acctrans_designations) > 1:
            i_o.split = True

        lineage_dict[intro] = i_o

    return lineage_objects, lineage_dict, unclear_taxa

def find_merged(lineage_objects):

    acctran_dict = defaultdict(list)
    acctran_to_merge = defaultdict(list)
    merged = {}
    merge_count = 0
    skip_list = []

    for lin in lineage_objects:         
        for i in lin.acctrans_designations: #should be only one
            acctran_dict[i].append(lin)
           
    for key, value in acctran_dict.items():
        if key not in skip_list:
            if len(value) > 1:
                for i in value:
                    merged[i.id] = True
                    merge_count += 1
                    if i not in acctran_to_merge[key]:
                        acctran_to_merge[key].append(i)

                    if i.split:
                        others = []
                        for j in i.acctrans_designations:
                            skip_list.append(j)
                            for k in acctran_dict[j]:
                                if k not in others and k not in acctran_to_merge[key]:
                                    others.append(k) #get the other lineages that need to go in
                        
                        acctran_to_merge[key].extend(others)
            else:
                if not i.split:
                    for i in value: #there will be one
                        merged[i.id] = False
        else:
            pass

    return merged, acctran_dict, acctran_to_merge, merge_count

def deal_with_merged(acctran_to_merge, lineage_object_dict, unclear_taxa, lineage_objects):

    remove_list = []
    keep_ids = []

    for acctran, lineages in acctran_to_merge.items():

        taxa = []
        numbers = []
        for_df = []
    
        for lin in lineages:
            taxa.extend(lin.taxa)
            remove_list.append(lin)
            for_df.append(lin.id)
            numbers.append(int(lin.id.lstrip("UK")))

        test_taxa = set(taxa)
        if len(test_taxa) != len(taxa):
            print("error in merging" + str(lineages) + " " + str(len(test_taxa)) + " " + str(len(taxa)))
        

        new_uk = "UK" + str(min(numbers))
        keep_ids.append(new_uk)

        for tax in unclear_taxa:
            if tax.acctrans == acctran:
                tax.lineage= new_uk
                taxa.append(tax)

        i_o = lineage(new_uk, taxa)
        i_o.acctrans_designations = acctran

        for tax in i_o.taxa:
            tax.lineage = i_o.id

        lineage_object_dict[i_o.id] = i_o      

        lineage_objects.append(i_o)

    for i in remove_list:
        if i in lineage_objects:
            lineage_objects.remove(i)

       
    return lineage_object_dict, lineage_objects


def write_to_file(lineage_objects, outfile):

    fw = open(outfile, 'w')
    fw.write("taxon,uk_lineage\n")
    for lin in lineage_objects:
        for tax in lin.taxa:
            fw.write(tax.id + "," + tax.lineage + "\n")

    fw.close()


def merge_lineages(input_dir, outfile):

    introduction_int_list, taxon_list, intros_to_taxa, acctrans_to_intro, intro_acctrans = make_taxon_objects(input_dir)

    lineage_objects, lineage_dict, unclear_taxa = deal_with_new_lineages(introduction_int_list, taxon_list, intros_to_taxa, acctrans_to_intro, intro_acctrans)

    merged, acctran_dict, acctran_to_merge, merge_count = find_merged(lineage_objects)

    lineage_object_dict, lineage_objects = deal_with_merged(acctran_to_merge, lineage_dict, unclear_taxa, lineage_objects)

    write_to_file(lineage_objects, outfile)