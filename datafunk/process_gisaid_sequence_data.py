from Bio import SeqIO
import json, re


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
    country_strings = gisaid_json_dict['covv_location'].split(" / ")[1:4]
    while len(country_strings) < 3:
        country_strings.append("")

    myStr = gisaid_json_dict['covv_virus_name'] + '|' + \
            gisaid_json_dict['covv_accession_id'] + '||' + \
            '|'.join(country_strings) + '|' + \
            gisaid_json_dict['covv_collection_date']

    return(myStr)

def update_fasta_header_string(header):
    fields = header.split('|')
    updated_fields = fields[:2]
    updated_fields.extend(["","UK"])
    updated_fields.extend(fields[2:])
    updated_header = '|'.join(updated_fields)
    return updated_header

def fix_header(header):
    """
    parse fasta header and remove problems
    """
    fixed_header = header.replace(' ', '_').replace("hCoV-19/","").replace("hCov-19/","")

    return(fixed_header)


def parse_omissions_file(file):
    """
    Parse a file of records to omit.

    Relies on there being a regex match to "EPI_ISL_\d{6}" in a line.

    Only returns the first match to the regex
    """
    # TO DO: allow regexes to e.g. 'bat', 'pangolin' in the omissions file

    IDs = []
    regex = re.compile('EPI_ISL_\d{6}')
    file_is_fasta = file.split('.')[-1][0:2] == 'fa'

    with open(file, 'r') as f:
        for line in f:
            if file_is_fasta:
                if line[0] != '>':
                    continue
            elif line.startswith("#"):
                continue
            match = re.search(regex, line)
            if match:
                ID = match.group()
                IDs.append(ID)

    return(IDs)

def keep_entry(header, omitted=False, exclude_uk=False):
    regex = re.compile('EPI_ISL_\d{6}')
    match = re.search(regex, header)
    if not match:
        return False
    epi_id = match.group()

    if omitted and epi_id in omitted:
        return False
    if exclude_uk:
        for country in ['/England/', '/Scotland/', '/Wales/', '/Northern Ireland/']:
            if country.lower() in header.lower():
                return False
    return True


def process_gisaid_sequence_data(input, output = False, omit_file_list = False, exclude_uk = False):

    def input_fasta_output(input, output, omitted, exclude_uk):
        if output:
            out = open(output, 'w')
        else:
            out = sys.stdout
        out = open(output, 'w')
        with open(input, 'r') as f:
            for record in SeqIO.parse(f, "fasta"):
                if keep_entry(record.description, omitted, exclude_uk):
                    out.write('>' + fix_header(update_fasta_header_string(record.description)) + '\n')
                    out.write(str(record.seq) + '\n')

        if output:
            out.close()
        pass

    def input_json_output(input, output, omitted, exclude_uk):
        if output:
            out = open(output, 'w')
        else:
            out = sys.stdout
        with open(input, 'r') as f:
                for jsonObj in f:
                    jsonDict = fix_seq_in_gisaid_json_dict(json.loads(jsonObj))
                    header = get_ID_from_json_dict(jsonDict)
                    if keep_entry(header, omitted, exclude_uk):
                        out.write('>' + fix_header(header) + '\n')
                        out.write(jsonDict['sequence'] + '\n')
        if output:
            out.close()
        pass


    if omit_file_list:
        temp = []
        for file in omit_file_list:
            temp.extend(parse_omissions_file(file))

        omitted_IDs = set(temp)
    else:
        omitted_IDs = False

    input_is_fasta = input.split('.')[-1][0:2].lower() == 'fa'
    input_is_json = input.split('.')[-1].lower() == 'json'
    if input_is_fasta:
        input_fasta_output(input = input, output = output, omitted = omitted_IDs, exclude_uk = exclude_uk)
    elif input_is_json:
        input_json_output(input = input, output = output, omitted = omitted_IDs, exclude_uk = exclude_uk)













#
