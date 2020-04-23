import pycountry
import pandas as pd
import os

df = pd.read_csv(os.path.dirname(os.path.realpath(__file__)) + '/resources/cities50000.tsv', sep = '\t')

cities_dict = {}
for i in range(len(df)):
    row = df.iloc[i,].squeeze()
    cities_dict[row.at['city'].lower()] = row.at['country_code']

countries_list = [x.name.lower() for x in pycountry.countries]

subdivisions_dict = {x.name.lower(): x.country_code for x in pycountry.subdivisions}

others = {'hubei': ('China', 'Hubei'),
          'wuhan': ('China', 'Hubei'),
          'korea': ('South_Korea', ''),
          'us': ('USA', ''),
          'u.s.a': ('USA', ''),
          'usa': ('USA', ''),
          'u.k': ('UK', ''),
          'uk': ('UK', ''),
          'gran.canary': ('Spain', 'Canary_Islands'),
          'irn': ('Iran', ''),
          'ny': ('USA', 'New_York_City'),
          'la': ('USA', 'Los_Angeles'),
          'iran': ('Iran', ''),
          'fareo': ('Faroe_Islands', ''),
          'czech': ('Czech_Republic', ''),
          'finnland': ('Finland', ''),
          'prague': ('Czech_Republic', 'Prague')}


def get_travel_history(json_dict):
    """
    travel history, if it exists, is in: json_dict['covv_add_host_info']
    """
    country = []
    subdivision = []
    city = []
    other = []

    add_host_info = json_dict['covv_add_host_info']
    add_host_info_split = add_host_info.split()

    for word in add_host_info_split:
        word = word.rstrip(",.:;!?-'\"").lstrip(",.:;!?-'\"")

        # These are valid cities, which will be excluded because they match other words:
        if word.lower() in ['of', 'holiday', 'asia', 'northern', 'sur']:
            continue

        if word.lower() in countries_list:
            country.append((word.replace(' ', '_'), ''))

        if word.lower() in subdivisions_dict:
            if word.lower() not in ['hubei', 'wuhan']:
                subdivision.append(word)
                sd_country = pycountry.countries.get(alpha_2=subdivisions_dict[word.lower()]).name
                if sd_country == 'Iran, Islamic Republic of':
                    country.append(('Iran', word))
                else:
                    country.append((sd_country.replace(' ', '_'), word))

        if word.lower() in cities_dict:
            if word.lower() not in ['hubei', 'wuhan', 'prague']:
                city.append(word)
                city_country = pycountry.countries.get(alpha_2=cities_dict[word.lower()]).name
                if city_country == 'Iran, Islamic Republic of':
                    country.append(('Iran', word))
                else:
                    country.append((city_country.replace(' ', '_'), word))

        if word.lower() in others:
            country.append(others[word.lower()])

    locales = set(country)

    # sometimes sampling country is represented in the text. Let's remove whereever
    # this is from the set, if it's in there
    sampling_country = json_dict['edin_admin_0'].replace(' ', '_')

    def test_sampling_country(item):
        return item[0] != sampling_country

    locales = set(filter(test_sampling_country, locales))

    country_only = [x for x in locales if len(x[1]) == 0]
    country_and_subdivision = [x for x in locales if len(x[1]) > 0]

    for country in country_only:
        if any([country[0] == x[0] for x in country_and_subdivision]):
            locales.remove(country)

    locales = [x[0] + '/' + x[1] if len(x[1]) > 0 else x[0] for x in locales]

    json_dict['edin_travel'] = ';'.join(locales).replace(',', '')

    return(json_dict)
