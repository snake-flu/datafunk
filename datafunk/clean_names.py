import pandas as pd
import pycountry as pc
import re

def clean_name(input_file, input_trait, output_file = "cleaned_file.csv"):
    log_file = open(output_file+".log","w")
    metadata = pd.read_csv(input_file)
    metadata.columns = metadata.columns.str.lower()
    trait = input_trait.lower()

    for index, row in metadata.iterrows():
        country = metadata.loc[index,trait]
        regex = re.compile('[^a-zA-Z]')
        search_item = regex.sub(' ',country.lower())
        try:
            search_result = pc.countries.search_fuzzy(search_item)
        except:
            log_file.write("Fail to parse country: " + str(country) + " for label ID: " + str(metadata.iloc[[index]]) + "\n")
        metadata.loc[index,trait] = search_result[0].name

    metadata.to_csv(output_file)    