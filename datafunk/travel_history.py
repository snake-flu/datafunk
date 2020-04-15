import pycountry
import pandas as pd

url = 'https://raw.githubusercontent.com/benjamincjackson/world_cities/master/cities50000.tsv'
df = pd.read_csv(url, sep = '\t')

cities_dict = {}
for i in range(len(df)):
    row = df.iloc[i,].squeeze()
    cities_dict[row.at['city'].lower()] = row.at['country_code']

countries_list = [x.name.lower() for x in pycountry.countries]

subdivisions_dict = {x.name.lower(): x.country_code for x in pycountry.subdivisions}

others = {'hubei': 'China',
          'wuhan': 'China',
          'korea': 'South Korea',
          'us': 'USA',
          'u.s.a': 'USA',
          'usa': 'USA',
          'u.k': 'UK',
          'uk': 'UK',
          'gran.canary': 'Spain',
          'irn': 'Iran',
          'ny': 'USA',
          'la': 'USA',
          'iran': 'Iran'}


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

        # Note these are valid cities, which will be excluded as it stands:
        if word.lower() in ['of', 'holiday', 'asia', 'northern']:
            continue

        if word.lower() in countries_list:
            country.append(word)
            continue

        if word.lower() in subdivisions_dict:
            subdivision.append(word)
            sd_country = pycountry.countries.get(alpha_2=subdivisions_dict[word.lower()]).name
            if sd_country == 'Iran, Islamic Republic of':
                country.append('Iran')
            else:
                country.append(sd_country)
            continue

        if word.lower() in cities_dict:
            city.append(word)
            city_country = pycountry.countries.get(alpha_2=cities_dict[word.lower()]).name
            country.append(city_country)
            continue

        if word.lower() in others:
            country.append(others[word.lower()])

    everything_else = set(subdivision + city + other)
    countries = set(country)

    # sometimes sampling country is represented in the text. Let's remove whereever
    # this is from the set, if it's in there
    sampling_country = json_dict['edin_admin_0']

    countries = countries - set([sampling_country])

    json_dict['edin_travel'] = ':'.join(countries)

    return(json_dict)
