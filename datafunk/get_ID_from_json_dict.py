import json

def get_ID_from_json_dict(gisaid_json_dict):
    """
    make a sequence identifier from the gisaid dump json-format data.

    Function input (gisaid_json_dict) is one record from the dump,
    e.g. do something like this:

    import json
    
    gisaid_dump = []
    with open('gisaid_dump.json', 'rt') as f:
        for jsonObj in f:
            gisaid_dump.append(json.loads(jsonObj))

    IDs = [get_ID_from_json_dict(x) for x in gisaid_dump]

    (note: 'Cov'/'CoV' both occur in the covv_virus_name field. These could be unified)

    """
    myStr = gisaid_json_dict['covv_virus_name'].replace(' ', '_') + '|' + \
            gisaid_json_dict['covv_accession_id'] + '|' + \
            gisaid_json_dict['covv_collection_date']
            
    return(myStr)
