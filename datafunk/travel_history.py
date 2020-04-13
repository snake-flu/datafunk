import pycountry
from datapackage import Package

"""
TO DO: lookup dict for admin names to unified output DONE
TO DO: exclude common words that are also city names ('of') PART
TO DO: include excluded cities because English spelling ('Milan') PART
TO DO: include common regions because not country, region, or city ('Hubei' = China) PART
TO DO: if 'Traveled from' but no match, flag for follow up
TO DO: maybe strip all punctuation from word fields? DONE
TO DO: other flags/log to output
TO DO: final format: ':'.join(set(places)) ? DONE
"""


package = Package('https://datahub.io/core/world-cities/datapackage.json')
c = package.get_resource('world-cities_csv')
cities_list = [x[0].lower() for x in c.read()]


extras_cities = ['ny', 'la', 'milan']
cities_list = cities_list + extras_cities


countries_list = [x.__getattr__('name').lower() for x in pycountry.countries]
extras_countries = ['iran', 'south korea', 'russia', 'korea', 'uk', 'usa', 'us', \
                    'u.k.', 'u.s.a.', 'gran.canary', 'irn']

countries_list = countries_list + extras_countries


subdivisions_list = [x.__getattr__('name').lower() for x in pycountry.subdivisions]


others = ['hubei', 'wuhan']


lookup = {'iran': 'Iran',
          'south korea': 'South Korea',
          'russia': 'Russia',
          'korea': 'Korea',
          'uk': 'UK',
          'usa': 'USA',
          'us': 'USA',
          'u.k': 'UK',
          'u.s.a': 'USA',
          'gran.canary': 'Gran Canaria',
          'irn': 'Iran',
          'ny': 'New York',
          'la': 'Los Angeles'}


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

        if word.lower() in lookup:
            word = lookup[word.lower()]
            # print(word)
            # if word.lower() in countries_list:
            #     country.append(word)
            #     print(country)
            #     print()
        # Note these are valid cities, which will be excluded as it stands:
        if word.lower() in ['of', 'holiday', 'asia', 'northern']:
            continue
        if word.lower() in countries_list:
            country.append(word)
        if word.lower() in subdivisions_list:
            subdivision.append(word)
        if word.lower() in cities_list:
            city.append(word)
        if word.lower() in others:
            other.append(word)

    everything = set(country + subdivision + city + other)

    json_dict['edin_travel'] = ':'.join(everything)

    return(json_dict)
