import requests
from astroquery.simbad import Simbad
import numpy as np
import pandas as pd
from astropy.table import QTable, Table, Column
from astropy import units as u
import urllib
import re
import bs4
import math
import matplotlib.pyplot as plt

def plot(star_name):

    # Convert the names of stars to HIP numbers, and get HIP numbers.
    if star_name.startswith('HIP'):
        HIP = int(star_name.split('HIP')[1])
    else:
        result_table = Simbad.query_objectids(star_name)
        line = list(filter(lambda x: 'HIP' in str(x), result_table))
        HIP = int(line[0][0].split('HIP')[1])

    try:
        url = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?hipepId={HIP}'
        webpage = str(urllib.request.urlopen(url).read())
        soup = bs4.BeautifulSoup(webpage,'html.parser')
        text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\r\\n'")
        text = text.split('\\n',17)
        label = [text[6],text[8],text[9]]
        label_list = [x.split(':',1)[1].strip().split('        ')[0] for x in label]
        label_list = [float(x) for x in label_list]
        print('Median Magnitude(red line) (mag):',label_list[0])
        print('5th percentile (max) (mag):',label_list[1])
        print('95th percentile (min) (mag):',label_list[2])
        text_list = text[17].split('\\r\\n')
        data_list = [x.split('|') for x in text_list]
        data_list_t = list(map(list, zip(*data_list)))
        index_float = [0,1,2]
        for i in index_float:
            data_list_t[i] = [float(x) for x in list(data_list_t[i])]

        text_list = text[17].split('\\r\\n')
        data_list = [x.split('|') for x in text_list]
        data_list_t = list(map(list, zip(*data_list)))
        index_float = [0,1,2]
        for i in index_float:
            data_list_t[i] = [float(x) for x in list(data_list_t[i])]

        plt.title(f'HIP {HIP}')
        plt.xlabel('Observation epoch (BJD - 2440000)')
        plt.ylabel('Calibrated Hp (mag)')
        plt.gca().invert_yaxis()
        plt.axhline(y=label_list[0],color='red')
        plt.scatter(data_list_t[0],data_list_t[1],s=10,color='blue')
        plt.errorbar(data_list_t[0],data_list_t[1],yerr=data_list_t[2],fmt='none',ecolor='cyan')

    except IndexError:
        print(f'The data of HIP {HIP} cannot be found.')