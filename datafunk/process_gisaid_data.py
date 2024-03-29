"""
write metadata and output sequences at the same time
"""

from Bio import SeqIO
import datetime
from datetime import datetime
from epiweeks import Week, Year
import sys
import json
import argparse
import warnings
import re
import pycountry
from itertools import chain
from unidecode import unidecode
import unicodedata

from datafunk.travel_history import *
from datafunk.travel_history import cities_dict, countries_list, subdivisions_dict, others

"""Don't edit these two lists please:
"""
_fields_gisaid = ['covv_accession_id', 'covv_virus_name', 'covv_location', 'covv_collection_date',
                 'covv_add_host_info', 'covv_assembly_method', 'covv_gender', 'covv_host',
                 'covv_passage', 'covv_patient_age', 'covv_seq_technology',
                 'covv_specimen', 'covv_subm_date',
                 'covv_patient_status', 'covv_lineage', 'covv_add_location', 'covv_clade']

_fields_edin = ['edin_header', 'edin_admin_0', 'edin_admin_1', 'edin_admin_2',
                'edin_travel', 'edin_date_stamp', 'edin_omitted', 'edin_epi_week', 'edin_epi_day',
                'edin_flag', 'is_uk']

"""You can edit this list:
"""
fields = ['sample_date', 'safe_sample_date', 'is_cog_uk']


"""Functions
"""

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def fix_gisaid_json_dict(gisaid_json_dict):
    """
    Remove commas from fields inside json dict

    And remove non-ascii characters
    """

    newDict = {}
    for x,y in gisaid_json_dict.items():
        newDict[x] = unidecode(str(y).replace(',', ''))

    return(newDict)

def clean_string(text, strip_symbols = False):
    # strip accents
    text = ''.join((c for c in unicodedata.normalize('NFD', text) if unicodedata.category(c) != 'Mn'))
    # replace spaces with _
    text = text.replace(' ', '_')
    
    if strip_symbols:
        p = re.compile(r'[^a-zA-Z0-9_/-]')
        text = p.sub('', text)
    return text

def get_admin_levels_from_json_dict(gisaid_json_dict, warnings = True):
    """
    get location strings from the gisaid location field
    use pycountry /
    """
    location_strings = [x.strip() for x in gisaid_json_dict['covv_location'].split("/")]

    while len(location_strings) < 4:
        location_strings.append("")

    continent = location_strings[0]
    country = location_strings[1]
    subdivision = location_strings[2]
    subsubdivision = location_strings[3]

    if any([x == country for x in ['England', 'Northern Ireland', 'Scotland', 'Wales']]):
        country = 'United Kingdom'
        subdivision = location_strings[1]
        subsubdivision = location_strings[2]

    if any([x == country for x in ['Alaska']]):
        country = 'USA'
        subdivision = location_strings[1]
        subsubdivision = location_strings[2]

    # some check here that there's a match to a real country
    # using pycountries?
    # First, these countries are known exceptions:
    if warnings:
        if all([country != x for x in ['Iran', 'South Korea', 'Russia', 'Korea', 'Democratic Republic of the Congo']]):
            try:
                pycountry.countries.lookup(country)
            except LookupError:
                eprint('Check country flagged for ' + gisaid_json_dict['covv_accession_id'] + \
                              '  ("' + country + '")')

                if len(gisaid_json_dict['edin_flag']) == 0:
                    gisaid_json_dict['edin_flag'] = 'check_country'
                elif len(gisaid_json_dict['edin_flag']) > 0:
                    gisaid_json_dict['edin_flag'] = gisaid_json_dict['edin_flag'] + ':check_country'

    if country == 'United Kingdom':
        country = 'UK'

    if country == 'Korea':
        country = 'South Korea'

    if country == 'Democratic Republic of the Congo':
        country = 'DRC'

    gisaid_json_dict['edin_admin_0'] = clean_string(country)
    gisaid_json_dict['edin_admin_1'] = clean_string(subdivision)
    gisaid_json_dict['edin_admin_2'] = clean_string(subsubdivision)

    return(gisaid_json_dict)


def expand_dict(dict, fields_list_required, fields_list_optional):
    """
    check fields in a dict
    make/delete them as appropriate
    """

    for x in fields_list_required:
        if x not in dict:
            dict[x] = ''

    for x in fields_list_optional:
        if x not in dict:
            dict[x] = ''

    return(dict)


def add_edin_flag_if_dict_key_true(dictionary, key):
    if key in dictionary:
        if dictionary[key] == "True" or dictionary[key] == "true" or dictionary[key] == True:
            if len(dictionary['edin_flag']) == 0:
                dictionary['edin_flag'] = key
            elif len(dictionary['edin_flag']) > 0 and key not in dictionary['edin_flag']:
                dictionary['edin_flag'] = dictionary['edin_flag'] + ':' + key
        del dictionary[key]
    return dictionary


def get_csv_order_and_record_dict(csv_file, fields_list_required, fields_list_optional):
    """
    Read all info in a csv metadata file into memory.

    old_records is a nested dict with EPI_IDs as the
    top-level keys and a dict of key: value pairs
    as the top-level values
    """
    first = True
    old_records = {}
    record_order = []
    with open(csv_file, 'r') as f:
        for line in f:
            l = line.strip().split(',')
            if first:
                keys = l
                extra_fields = list(set(keys) - set(fields_list_required + fields_list_optional))

                first = False
                continue

            if len(l) != len(keys):
                sys.exit('Badly formatted csv file: ' + csv_file.csv + ', exiting program')

            d = {x: y for x,y in zip(keys, l)}

            ID = d['covv_accession_id']
            record_order.append(ID)
            d = add_edin_flag_if_dict_key_true(d, "subsample_omit")
            old_records[ID] = d

    return(record_order, old_records, extra_fields)


def get_json_order_and_record_dict(json_file, fields_list_required, fields_list_optional):
    """
    Read all info in a GISAID json dump into memory.

    all_records is a nested dict with EPI_IDs as the
    top-level keys and a dict of key: value pairs
    as the top-level values
    """
    all_records = {}
    record_order = []
    extra_fields = []
    with open(json_file, 'r') as f:
        for jsonObj in f:

            d = fix_gisaid_json_dict(json.loads(jsonObj))
            d = fix_seq_in_gisaid_json_dict(d)

            extra_fields = extra_fields + list(set(list(d.keys())) - set(fields_list_required + fields_list_optional + extra_fields))

            ID = d['covv_accession_id']
            record_order.append(ID)
            all_records[ID] = d

    return(record_order, all_records, extra_fields)


def check_gisaid_date(dict):
    """
    Exclude dates that are incomplete or impossible (in the future)
    """
    date = dict['covv_collection_date']
    regex = re.compile('\d{4}-\d{2}-\d{2}')
    match = re.search(regex, date)
    if not match:
        # dict['edin_omitted'] = 'True'

        if len(dict['edin_flag']) == 0:
            dict['edin_flag'] = 'omitted_date'
        elif len(dict['edin_flag']) > 0:
            dict['edin_flag'] = dict['edin_flag'] + ':omitted_date'

    else:
        today = datetime.today().date()
        gisaid_date_as_python_date = datetime.strptime(date, '%Y-%m-%d').date()
        if gisaid_date_as_python_date > today or gisaid_date_as_python_date < datetime.strptime('2019-11-30','%Y-%m-%d').date():
            if len(dict['edin_flag']) == 0:
                dict['edin_flag'] = 'omitted_date'
            elif len(dict['edin_flag']) > 0:
                dict['edin_flag'] = dict['edin_flag'] + ':omitted_date'

    return(dict)


def check_edin_omitted_file(dict, omit_set):
    """
    update omission field to True if record
    is in omissions file
    """
    if omit_set:
        if dict['covv_accession_id'] in omit_set:
            # dict['edin_omitted'] = 'True'

            if len(dict['edin_flag']) == 0:
                dict['edin_flag'] = 'omitted_file'
            elif len(dict['edin_flag']) > 0:
                dict['edin_flag'] = dict['edin_flag'] + ':omitted_file'

        elif any([x in dict['covv_virus_name'] for x in ['/bat/', '/pangolin/']]):

            if len(dict['edin_flag']) == 0:
                dict['edin_flag'] = 'omitted_file'
            elif len(dict['edin_flag']) > 0:
                dict['edin_flag'] = dict['edin_flag'] + ':omitted_file'

        elif dict['covv_host'].lower() != 'human':

            if len(dict['edin_flag']) == 0:
                dict['edin_flag'] = 'omitted_file'
            elif len(dict['edin_flag']) > 0:
                dict['edin_flag'] = dict['edin_flag'] + ':omitted_file'

        elif any([x in dict['edin_header'] for x in omit_set]):

            if len(dict['edin_flag']) == 0:
                dict['edin_flag'] = 'omitted_file'
            elif len(dict['edin_flag']) > 0:
                dict['edin_flag'] = dict['edin_flag'] + ':omitted_file'



    return(dict)


def update_UK_sequence(gisaid_json_dict):
    """
    Flag UK sequences
    """
    header = gisaid_json_dict['edin_header']
    gisaid_json_dict["is_cog_uk"] = False
    gisaid_json_dict["is_uk"] = False
    for country in ['England/', 'Scotland/', 'Wales/', 'Northern_Ireland/']:
        if country.lower() in header.lower():

            if len(gisaid_json_dict['edin_flag']) == 0:
                gisaid_json_dict['edin_flag'] = 'uk_sequence'
            elif len(gisaid_json_dict['edin_flag']) > 0:
                gisaid_json_dict['edin_flag'] = gisaid_json_dict['edin_flag'] + ':uk_sequence'

            gisaid_json_dict['is_uk'] = 'True'
            gisaid_json_dict['is_cog_uk'] = 'True'

    for country in ['Falkland_Islands/', 'Gibraltar/', 'Jersey/', 'Isle_of_Man/', 'Guernsey/']:
        if country.lower() in header.lower():
            gisaid_json_dict['is_cog_uk'] = 'True'
    return(gisaid_json_dict)


def update_edin_date_stamp_field(gisaid_json_dict):
    """
    date stamp the new records to write
    """
    mydate = str(datetime.date(datetime.now()))
    gisaid_json_dict['edin_date_stamp'] = mydate
    return(gisaid_json_dict)


def date_string_to_epi_week(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    cumulative epi week which is cumulative total epidemiological
    weeks since 2019-12-22. Week beginning 2019-12-22 is week 0
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week:
    week = Week.fromdate(date)
    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    elif week.year == 2019:
        return("0")
    else:
        cum_epi_week = week.week + len(list(chain(*[[x for x in Year(y).iterweeks()] for y in range(2020, week.year)])))
        return str(cum_epi_week)


def date_string_to_epi_day(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    cumulative epi day which is cumulative total days since 2019-12-22
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week week:
    week = Week.fromdate(date)
    # this is day 1 of epi-week 0:
    day_one = datetime.strptime("2019-12-22", '%Y-%m-%d').date()
    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    else:
        cum_epi_day = (date - day_one).days + 1
        return str(cum_epi_day)

def date_string_to_safe_date_string(date_string):
    """
    parse a date string in YYYY-MM-DD format and return
    date corresponding to the start of the epi-week in which it falls.
    Week beginning 2019-12-22 is week 0
    """
    try:
        date = datetime.strptime(date_string, '%Y-%m-%d').date()
    except:
        return ""
    # this is epi-week:
    week = Week.fromdate(date)

    if week.year < 2019 or (week.year == 2019 and week.week < 52):
        return ""
    else:
        return week.startdate().strftime('%Y-%m-%d')

def update_edin_epi_date_fields(gisaid_json_dict):
    """
    record epi week and epi day by parsing sample
    collection date

    NB this will break in 2021!
    """
    if 'omitted_date' in gisaid_json_dict['edin_flag']:
        return(gisaid_json_dict)

    collection_date = gisaid_json_dict['covv_collection_date']

    # returns None if nothing found
    epi_week = date_string_to_epi_week(collection_date)
    epi_day = date_string_to_epi_day(collection_date)
    safe_date = date_string_to_safe_date_string(collection_date)

    if epi_week:
        gisaid_json_dict['edin_epi_week'] = epi_week

    if epi_day:
        gisaid_json_dict['edin_epi_day'] = epi_day

    if safe_date:
        gisaid_json_dict['sample_date'] = collection_date
        gisaid_json_dict['safe_sample_date'] = safe_date

    return(gisaid_json_dict)


def get_one_metadata_line(dict, fields_list, sep = ','):
    """
    return a string of dict.values() formatted
    for writing in the order in fields_list

    l is a 'sep'-separated string formatted for printing
    """
    l = sep.join([str(dict[x]) for x in fields_list]) + '\n'
    return(l)


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


def fix_header(header):
    """
    parse fasta header and remove problems
    """
    fixed_header = header.lstrip()\
        .replace(' ', '_')\
        .replace("hCoV-19/","")\
        .replace("hCov-19/","")\
        .replace("PENDING", "")\
        .replace("UPLOADED", "")\
        .replace("None", "")

    fixed_header = re.sub('^Wuhan-Hu-1','China/Wuhan-Hu-1',fixed_header)
    fixed_header = clean_string(fixed_header)

    return(fixed_header)


def add_header_to_json_dict(gisaid_json_dict):
    """
    make a sequence identifier from the gisaid dump json-format data.
    Function input (gisaid_json_dict) is one record from the dump
    """

    myTempHeader = gisaid_json_dict['covv_virus_name'] + '|' + \
                    gisaid_json_dict['covv_accession_id'] + '||' + \
                    gisaid_json_dict['edin_admin_0'] + '|' + \
                    gisaid_json_dict['edin_admin_1'] + '|' + \
                    gisaid_json_dict['edin_admin_2'] + '|' + \
                    gisaid_json_dict['covv_collection_date']

    gisaid_json_dict['edin_header'] = fix_header(myTempHeader)

    return(gisaid_json_dict)


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
            elif len(line.strip()) == 0:
                continue

            match = re.search(regex, line)
            if match:
                ID = match.group()
                IDs.append(ID)
            else:
                IDs.append(line.rstrip())

    return(IDs)


def write_metadata_output(output,
                  new_records_list,
                  new_records_dict,
                  old_records_list,
                  old_records_dict,
                  fields_list):

    """
    write a csv-format outfile to file
    """
    out = open(output, 'w')

    out.write(','.join(fields_list) + '\n')

    for record in old_records_list:
        do = old_records_dict[record]
        lo = get_one_metadata_line(dict=do, fields_list=fields_list)
        out.write(lo)

    for record in new_records_list:
        dn = new_records_dict[record]
        ln = get_one_metadata_line(dict=dn, fields_list=fields_list)
        out.write(ln)

    out.close()
    pass


def write_fasta_output(output,
                      new_records_list,
                      new_records_dict,
                      old_records_list,
                      old_records_dict):
    """
    write the sequences to a fasta file
    """
    if output:
        out = open(output, 'w')
    else:
        out = sys.stdout

    for record in old_records_list:
        if old_records_dict[record]['edin_omitted'] == 'True':
            continue

        else:
            header = old_records_dict[record]['edin_header']
            seq = old_records_dict[record]['sequence']
            out.write('>' + header + '\n')
            out.write(seq + '\n')


    for record in new_records_list:
        if new_records_dict[record]['edin_omitted'] == 'True':
            continue

        else:
            header = new_records_dict[record]['edin_header']
            seq = new_records_dict[record]['sequence']
            out.write('>' + header + '\n')
            out.write(seq + '\n')


    if output:
        out.close()
    pass


def compare_records(metadata_gisaid_dict, json_gisaid_dict):
    """
    check for equality between gisaid fields in the new dump
    vs. the last iteration of the metadata. If equality == False
    then we reprocess this record
    """
    equality = True
    for field in _fields_gisaid:
        a = metadata_gisaid_dict[field]
        b = json_gisaid_dict[field]
        if a != b:
            equality = False
    return(equality)


def repopulate_sequence_from_new_dump(csv_record_dict, all_records_dict):
    """
    add sequences back into old references from the new dump
    """
    # this is the id of the old record
    epi_id = csv_record_dict['covv_accession_id']
    # this is the sequence for this id in the new dump
    seq = all_records_dict[epi_id]['sequence']
    # add and populate the sequence field
    csv_record_dict['sequence'] = seq
    # return it
    return(csv_record_dict)


def wipe_edin_omit_field(json_gisaid_dict):
    json_gisaid_dict['edin_omitted'] = ''
    return(json_gisaid_dict)


def update_edin_omit_field(json_gisaid_dict,
                           exclude_uk = False, exclude_undated = False,
                           exclude_subsampled = True, exclude_omitted_file = True):

    if exclude_uk:
        if 'uk_sequence' in json_gisaid_dict['edin_flag']:
            json_gisaid_dict['edin_omitted'] = 'True'
            return(json_gisaid_dict)

    if exclude_undated:
        if 'omitted_date' in json_gisaid_dict['edin_flag']:
            json_gisaid_dict['edin_omitted'] = 'True'
            return(json_gisaid_dict)

    if exclude_subsampled:
        if 'subsample_omit' in json_gisaid_dict['edin_flag']:
            json_gisaid_dict['edin_omitted'] = 'True'
            return(json_gisaid_dict)

    if exclude_omitted_file:
        if 'omitted_file' in json_gisaid_dict['edin_flag']:
            json_gisaid_dict['edin_omitted'] = 'True'
            return(json_gisaid_dict)

    return(json_gisaid_dict)


"""Program
"""
def process_gisaid_data(input_json,
                        input_omit_file_list,
                        input_metadata,
                        output_fasta,
                        output_metadata,
                        exclude_uk,
                        exclude_undated,
                        exclude_subsampled,
                        exclude_omitted_file):

    # logfile = open(output + '.log', 'w')
    if input_omit_file_list:
        temp = []
        for file in input_omit_file_list:
            temp.extend(parse_omissions_file(file))

        omitted_IDs = set(temp)
    else:
        omitted_IDs = False


    if input_metadata != 'False':
        # Check that all required fields were in the csv file
        csv_header = next(open(input_metadata, 'r')).strip().split(',')
        if not all([x in csv_header for x in _fields_gisaid + _fields_edin]):
            sys.exit('There were missing mandatory fields in ' + input_metadata)

        # # add optional fields from the old csv file to the output
        # for x in csv_header:
        #     if x not in _fields_gisaid + _fields_edin:
        #         if x not in fields:
        #             fields.append(x)


        old_records = get_csv_order_and_record_dict(input_metadata,
                                                    fields_list_required = _fields_gisaid + _fields_edin,
                                                    fields_list_optional = fields)

        temp_old_records_list = old_records[0]
        temp_old_records_dict = old_records[1]
        extra_fields = old_records[2]

        for x in extra_fields:
            fields.append(x)


    else:
        temp_old_records_list = []
        old_records_dict = {}


    all_records = get_json_order_and_record_dict(input_json,
                                                fields_list_required = _fields_gisaid + _fields_edin,
                                                fields_list_optional = fields)

    all_records_list = all_records[0]
    all_records_dict = all_records[1]
    new_fields = set(all_records[2]) - set(['sequence_length', 'sequence'])

    if len(new_fields) > 0:
        eprint('new fields in gisaid dump: ' + '\t'.join(new_fields))


    # for each old record:
    # if the info in the csv doesn't match the info in the new json dump,
    # throw this record out of the list of old records - which means that
    # it will get re-processed

    to_remove = []
    for record in temp_old_records_list:
        # this record might have been removed.
        # so check if it is in all_records_dict before proceeding
        if record not in all_records_dict:
            to_remove.append(record)
        elif temp_old_records_list.count(record) > 1:
            to_remove.append(record)
        else:
            old_record = temp_old_records_dict[record]
            new_record = all_records_dict[record]
            is_it_the_same = compare_records(old_record, new_record)
            if is_it_the_same == False:
                # this could go into a log:
                # print('\t'.join([old_record[x] for x in _fields_gisaid]))
                # print('\t'.join([new_record[x] for x in _fields_gisaid]))
                # print()
                to_remove.append(record)
    # this could go into a log:
    # print('removed because different: '+ str(changecount_diff))
    # print('removed because deleted: '+ str(changecount_del))
    # print('removed old = '+str(len(to_remove)))

    old_records_list = [x for x in temp_old_records_list if x not in to_remove]


    if input_metadata != 'False':

        # repopulate the old records with sequence from the new dump:
        old_records_dict = {x: repopulate_sequence_from_new_dump(temp_old_records_dict[x], all_records_dict) for x in set(old_records_list)}

        # expand dict to include any extra columns
        old_records_dict = {x: expand_dict(old_records_dict[x], fields_list_required = _fields_edin + _fields_gisaid, fields_list_optional = fields) for x in old_records_dict.keys()}

        # FIRST THING TO DO: WIPE EDIN_OMITTED
        old_records_dict = {x: wipe_edin_omit_field(old_records_dict[x]) for x in old_records_dict.keys()}

        # # TEMPORARY THINGS TO DO TO BRING OLD METADATA INLINE WITH NEW METADATA:
        # old_records_dict = {x: get_admin_levels_from_json_dict(old_records_dict[x]) for x in old_records_dict.keys()}
        # old_records_dict = {x: get_travel_history(old_records_dict[x]) for x in old_records_dict.keys()}

        # update omit field for this round of writing records only:
        old_records_dict = {x: update_edin_omit_field(old_records_dict[x],
                                           exclude_uk = exclude_uk,
                                           exclude_undated = exclude_undated,
                                           exclude_subsampled = exclude_subsampled,
                                           exclude_omitted_file = exclude_omitted_file)
                                for x in old_records_dict.keys()}



    # get new records out of the new dump:
    new_records_list = [x for x in all_records_list if x not in set(old_records_list)]
    new_records_dict = {x: all_records_dict[x] for x in new_records_list}

    # expand dict to include any extra columns
    new_records_dict = {x: expand_dict(new_records_dict[x], fields_list_required = _fields_edin + _fields_gisaid, fields_list_optional = fields) for x in new_records_dict.keys()}

    # FIRST THING TO DO: WIPE EDIN_OMITTED
    new_records_dict = {x: wipe_edin_omit_field(new_records_dict[x]) for x in new_records_dict.keys()}

    # update date stamp field
    new_records_dict = {x: update_edin_date_stamp_field(new_records_dict[x]) for x in new_records_dict.keys()}

    # update admin level
    new_records_dict = {x: get_admin_levels_from_json_dict(new_records_dict[x]) for x in new_records_dict.keys()}

    # include a header field in each dictionary (just for writing the fasta file):
    new_records_dict = {x: add_header_to_json_dict(new_records_dict[x]) for x in new_records_dict.keys()}

    # get travel history
    new_records_dict = {x: get_travel_history(new_records_dict[x]) for x in new_records_dict.keys()}

    # if gisaid collection date formatted correctly, we can add epi week and epi day
    new_records_dict = {x: update_edin_epi_date_fields(new_records_dict[x]) for x in new_records_dict.keys()}

    # check gisaid collection date formatted correctly
    new_records_dict = {x: check_gisaid_date(new_records_dict[x]) for x in new_records_dict.keys()}

    # check if sequence is in omissions file
    new_records_dict = {x: check_edin_omitted_file(new_records_dict[x], omitted_IDs) for x in new_records_dict.keys()}

    # record if sequence is from the UK
    new_records_dict = {x: update_UK_sequence(new_records_dict[x]) for x in new_records_dict.keys()}

    # update omit field for this round of writing records only:
    new_records_dict = {x:
        update_edin_omit_field(new_records_dict[x],
                                   exclude_uk = exclude_uk,
                                   exclude_undated = exclude_undated,
                                   exclude_subsampled = exclude_subsampled,
                                   exclude_omitted_file = exclude_omitted_file)
                        for x in new_records_dict.keys()}

    if output_metadata:
        write_metadata_output(output = output_metadata,
                              new_records_list = new_records_list,
                              new_records_dict = new_records_dict,
                              old_records_list = old_records_list,
                              old_records_dict = old_records_dict,
                              fields_list = _fields_edin + fields + _fields_gisaid)


    write_fasta_output(output = output_fasta,
                       new_records_list = new_records_list,
                       new_records_dict = new_records_dict,
                       old_records_list = old_records_list,
                       old_records_dict = old_records_dict)



    # logfile.close()
    pass




















#
