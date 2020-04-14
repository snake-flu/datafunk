import pycountry
from datapackage import Package

"""
TO DO: a better cities lookup dict (probably want major cities only)
"""


package = Package('https://datahub.io/core/world-cities/datapackage.json')
c = package.get_resource('world-cities_csv')
cities_dict = {x[0].lower(): {'city': x[0],
                              'country': x[1],
                              'subdivision': x[2],
                              'geoname': x[3]} for x in c.read()}

cities_dict['ny'] = {'city': 'New York',
                     'country': 'USA',
                     'subdivision': 'New York',
                     'geoname': 5128581}

cities_dict['la'] = {'city': 'Los Angeles',
                     'country': 'USA',
                     'subdivision': 'California',
                     'geoname': 5368361}

cities_dict['milan'] = {'city': 'Milan',
                     'country': 'Italy',
                     'subdivision': 'Lomardy',
                     'geoname': 3173435}

countries_list = [x.__getattr__('name').lower() for x in pycountry.countries]

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

        if word.lower() in subdivisions_dict:
            subdivision.append(word)
            sd_country = pycountry.countries.get(alpha_2=subdivisions_dict[word.lower()]).name
            if sd_country == 'Iran, Islamic Republic of':
                country.append('Iran')
            else:
                country.append(sd_country)

        if word.lower() in cities_dict:
            city.append(word)
            if word.lower() == 'venice' or word.lower() == 'rome':
                country.append('Italy')
            else:
                city_country = cities_dict[word.lower()]['country']
                country.append(city_country)

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
