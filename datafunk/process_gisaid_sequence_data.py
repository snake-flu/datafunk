from Bio import SeqIO
import json, re, argparse


def fix_seq_in_gisaid_json_dict(gisaid_json_dict):
    """
    strip whitespace and newline characters from the
    ['sequence'] field of a gisaid json object
    """

    def fix_seq(seq):
        newseq = re.sub(r"\s+", '', seq).replace('\n','')
        return(newseq)

    newDict = {}
    for x,y in gisaid_json_dict.items():
        if x == 'sequence':
            newDict['sequence'] = fix_seq(y)
        else:
            newDict[x] = y

    return(newDict)


def get_ID_from_json_dict(gisaid_json_dict):
    """
    make a sequence identifier from the gisaid dump json-format data.
    Function input (gisaid_json_dict) is one record from the dump
    """
    myStr = gisaid_json_dict['covv_virus_name'].replace(' ', '_') + '|' + \
            gisaid_json_dict['covv_accession_id'] + '|' + \
            gisaid_json_dict['covv_collection_date']

    return(myStr)


def parse_omissions_file(file):
    """
    Parse a file of records to omit.

    Relies on there being a regex match to "EPI_ISL_\d{6}" in a line.

    Only returns the first match to the regex
    """
    # TO DO: allow regexes to e.g. 'bat', 'pangolin' in the omissions file

    IDs = []
    regex = re.compile('EPI_ISL_\d{6}')

    with open(file, 'r') as f:
        for line in f:
            # if file is a FASTA file, only search the header lines:
            if file.split('.')[-1][0:2] == 'fa':
                if line[0] != '>':
                    continue
            match = re.search(regex, line)
            if match:
                ID = match.group()
                IDs.append(ID)

    return(IDs)


def process_gisaid_sequence_data(input, output = False, omit_file_list = False):

    def input_fasta_output_file(input, output, omitted):
        out = open(output, 'w')
        with open(input, 'r') as f:
            if omitted:
                omitted_IDs = omitted
                for record in SeqIO.parse(f, "fasta"):
                    # To do: try except IndexError below, in case header line not formatted correctly
                    ID = record.description.split('|')[1]
                    if ID not in omitted_IDs:
                        out.write('>' + record.description.replace(' ', '_') + '\n')
                        out.write(str(record.seq) + '\n')
            else:
                for record in SeqIO.parse(f, "fasta"):
                    out.write('>' + record.description.replace(' ', '_') + '\n')
                    out.write(str(record.seq) + '\n')
        out.close()
        pass


    def input_fasta_output_stdout(input, omitted):
        with open(input, 'r') as f:
            if omitted:
                omitted_IDs = omitted
                for record in SeqIO.parse(f, "fasta"):
                    # To do: try except IndexError below, in case header line not formatted correctly
                    ID = record.description.split('|')[1]
                    if ID not in omitted_IDs:
                        print('>' + record.description.replace(' ', '_'))
                        print(str(record.seq))
            else:
                for record in SeqIO.parse(f, "fasta"):
                    print('>' + record.description.replace(' ', '_'))
                    print(str(record.seq))
        pass


    def input_json_output_file(input, output, omitted):
        out = open(output, 'w')
        with open(input, 'r') as f:
            if omitted:
                for jsonObj in f:
                    jsonDict = fix_seq_in_gisaid_json_dict(json.loads(jsonObj))
                    ID = jsonDict['covv_accession_id']
                    if ID not in omitted_IDs:
                        out.write('>' + get_ID_from_json_dict(jsonDict) + '\n')
                        out.write(jsonDict['sequence'] + '\n')
            else:
                for jsonObj in f:
                    jsonDict = fix_seq_in_gisaid_json_dict(json.loads(jsonObj))
                    out.write('>' + get_ID_from_json_dict(jsonDict) + '\n')
                    out.write(jsonDict['sequence'] + '\n')
        out.close()
        pass


    def input_json_output_stdout(input, omitted):
        with open(input, 'r') as f:
            if omitted:
                for jsonObj in f:
                    jsonDict = fix_seq_in_gisaid_json_dict(json.loads(jsonObj))
                    ID = jsonDict['covv_accession_id']
                    if ID not in omitted_IDs:
                        print('>' + get_ID_from_json_dict(jsonDict))
                        print(jsonDict['sequence'])
            else:
                for jsonObj in f:
                    jsonDict = fix_seq_in_gisaid_json_dict(json.loads(jsonObj))
                    print('>' + get_ID_from_json_dict(jsonDict))
                    print(jsonDict['sequence'])
        pass


    if omit_file_list:
        temp = []
        for file in omit_file_list:
            temp.extend(parse_omissions_file(file))

        omitted_IDs = set(temp)
    else:
        omitted_IDs = None


    if input.split('.')[-1][0:2].lower() == 'fa':
        if output == 'stdout':
            input_fasta_output_stdout(input = input, omitted = omitted_IDs)
        else:
            if output:
                input_fasta_output_file(input = input, output = output, omitted = omitted_IDs)

    if input.split('.')[-1].lower() == 'json':
        if output == 'stdout':
            input_json_output_stdout(input = input, omitted = omitted_IDs)
        else:
            if output:
                input_json_output_file(input = input, output = output, omitted = omitted_IDs)

    pass













#
