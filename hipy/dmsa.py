import requests
from astroquery.simbad import Simbad
import numpy as np
import pandas as pd
from astropy.table import QTable, Table, Column
from astropy import units as u
import urllib
import re
import bs4

def get_dmsa(star_name):
    if star_name.startswith('HIP'):
        HIP = int(star_name.split('HIP')[1])
    else:
        result_table = Simbad.query_objectids(star_name)
        line = list(filter(lambda x: 'HIP' in str(x), result_table))
        HIP = int(line[0][0].split('HIP')[1])

    try:
        url1 = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?hipId={HIP}'
        webpage = str(urllib.request.urlopen(url1).read()) # Load the webpage
        # Parse webpage using https://www.crummy.com/software/BeautifulSoup/bs4/doc/#quick-start
        soup = bs4.BeautifulSoup(webpage, features='html.parser')
        # Get the text from the webpage. Note that the table is under the "pre" tag
        text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\n'") # remove extra '\\n' 
        text_list = text.split('\\n')
        flag = text_list[60].split(':')[1].lstrip(' ').split('                    ')[0]
    except IndexError:
        print(f'The catalogue of HIP {HIP} cannot be found.')
    
    url2 = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?dmId={HIP}'
    webpage = str(urllib.request.urlopen(url2).read())
    soup = bs4.BeautifulSoup(webpage, features='html.parser')
    text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\n'") 
    text_list = text.split('\\n')
    
    if flag == 'C':
        print(f'HIP {HIP} is in a component solution.')
        
        
    elif flag == 'G':
        print(f'HIP {HIP} is in an acceleration solution.')
        data_list = [x.split(':',1)[1].lstrip() for x in text_list[1:14]]
        value_list = [re.compile("\s{3,}").split(x) for x in data_list]
        void_flag = np.array([len(x) == 1 for x in value_list])
        _ = [value_list[i].insert(0, np.nan) for i in np.where(void_flag)[0]]
        value_list_t = list(zip(*value_list))
        final_list = [value_list_t[0]]
        corr_list = text_list[15] + text_list[16]
        final_list = np.append(final_list,[[corr_list]],1)
        final_list_t = list(zip(*final_list))
        final_list_t[1] = [float(x) for x in list(final_list_t[1])] * u.mas/(u.yr)**2
        final_list_t[2] = [float(x) for x in list(final_list_t[2])] * u.mas/(u.yr)**2
        final_list_t[3] = [float(x) for x in list(final_list_t[3])] * u.mas/(u.yr)**2
        final_list_t[4] = [float(x) for x in list(final_list_t[4])] * u.mas/(u.yr)**2
        final_list_t[6] = [float(x) for x in list(final_list_t[6])] * u.mas/(u.yr)**3
        final_list_t[7] = [float(x) for x in list(final_list_t[7])] * u.mas/(u.yr)**3
        final_list_t[8] = [float(x) for x in list(final_list_t[8])] * u.mas/(u.yr)**3
        final_list_t[9] = [float(x) for x in list(final_list_t[9])] * u.mas/(u.yr)**3
        out = QTable(final_list_t, names=('hip','g_ra*','g_dec','g_ra*_error','g_dec_error','f_g','gdot_ra*','gdot_dec',
                                  'gdot_ra*_error','gdot_dec_error','f_gdot','notes_flag','n','corr'),
            meta={'hip':'Hipparcos Catalogue identifier','g_ra*':'Components in RA of the apparent acceleration of the potocentre at epoch J1991.25 (mas/yr^2)',
                 'g_dec':'Components in Dec of the apparent acceleration of the potocentre at epoch J1991.25 (mas/yr^2)',
                 'g_ra*_error':'Standard error of g_ra* (mas/yr^2)','g_dec_error':'Standard error of g_dec (mas/yr^2)','f_g':'Significance of the g terms',
                 'gdot_ra*':'Components in RA of the rate of change of the apparent acceleration of the photocentre (mas/yr^3)',
                 'gdot_dec':'Components in RA of the rate of change of the apparent acceleration of the photocentre (mas/yr^3)',
                 'gdot_ra*_error':'Standard error of g_dot_ra* (mas/yr^3)','gdot_dec_error':'Standard error of g_dot_dec (mas/yr^3)',
                  'f_gdot':'Significance of the gdot terms','notes_flag':'Flag indicating a note at the end of relevent volume',
                 'n':'Number of astrometric parameters, n = 7 or 9','corr':'(Up to) 36 correlation coefficients in a non-linear coding'})

    elif flag == 'O':
        print(f'HIP {HIP} is in an orbital solution.')
        data_list = [x.split(':',1)[1].lstrip() for x in text_list[1:19]]
        value_list = [re.compile("\s{3,}").split(x) for x in data_list]
        void_flag = np.array([len(x) == 1 for x in value_list])
        _ = [value_list[i].insert(0, np.nan) for i in np.where(void_flag)[0]]
        value_list_t = list(zip(*value_list))
        final_list = [value_list_t[0]]
        corr_list = text_list[20] + text_list[21] + text_list[22] + text_list[23]
        final_list = np.append(final_list,[[corr_list]],1)
        final_list_t = list(zip(*final_list))
        final_list_t[1] = [float(x) for x in list(final_list_t[1])] * u.d
        final_list_t[2] = [float(x)+2440000 for x in list(final_list_t[2])] * u.d
        final_list_t[3] = [float(x) for x in list(final_list_t[3])] * u.mas
        final_list_t[4] = [float(x) for x in list(final_list_t[4])]
        final_list_t[5] = [float(x) for x in list(final_list_t[5])] * u.deg
        final_list_t[6] = [float(x) for x in list(final_list_t[6])] * u.deg
        final_list_t[7] = [float(x) for x in list(final_list_t[7])] * u.deg
        final_list_t[8] = [float(x) for x in list(final_list_t[8])] * u.d
        final_list_t[9] = [float(x) for x in list(final_list_t[9])] * u.d
        final_list_t[10] = [float(x) for x in list(final_list_t[10])] * u.mas
        final_list_t[11] = [float(x) for x in list(final_list_t[11])]
        final_list_t[12] = [float(x) for x in list(final_list_t[12])] * u.deg
        final_list_t[13] = [float(x) for x in list(final_list_t[13])] * u.deg
        final_list_t[14] = [float(x) for x in list(final_list_t[14])] * u.deg
        out = QTable(final_list_t, names=('hip','p','t','a0','e','periastron_argument','i','node_position_angle',
                                  'p_error','t_error','a0_error','e_error','periastron_argument_error','i_error','node_position_angle_error',
                                 'ref','notes_flag','status_flags','corr'),
            meta={'hip':'Hipparcos Catalogue identifier','p':'Orbital period, in days','t':'Time of periastron passage, in days',
                  'a0':'Semi-major axis of photocentre orbit (mas)','e':'Eccentricity','periastron_argument':'Argument of periastron (deg)',
                  'i':'Inclination (deg)','node_position_angle':'Position angle of the node (deg)',
                  'p_error':'Standard of error of p (days)','t_error':'Standard error of t (days)','a0_error':'Standard error of a0 (mas)',
                  'e_error':'Standard error of e','periastron_argument_error':'Standard error of the argument of periastron',
                  'i_error':'Standard error of i','node_position_angle_error':'Standard error of the position angle of the node',
                  'ref':'Reference to the literature for the orbital parameters','notes_flag':'Flag indicating a note at the end of relevent volume',
                  'status_flags':'Status flags for the 12 astrometric and orbital parameters taken in the order indicated below (1 = estimated, 0 = not estimated)',
                  'corr':'66 correlation coefficients in a non-linear coding'})

    elif flag == 'V':
        print(f'HIP {HIP} is in a solution of variability-induced movers.')
        data_list = [x.split(':',1)[1].lstrip() for x in text_list[1:13]]
        value_list = [re.compile("\s{3,}").split(x) for x in data_list]
        void_flag = np.array([len(x) == 1 for x in value_list])
        _ = [value_list[i].insert(0, np.nan) for i in np.where(void_flag)[0]]
        value_list_t = list(zip(*value_list))
        final_list = [value_list_t[0]]
        corr_list = text_list[14] + text_list[15]
        final_list = np.append(final_list,[[corr_list]],1)
        final_list_t = list(zip(*final_list))
        final_list_t[1] = [float(x) for x in list(final_list_t[1])] * u.mag
        final_list_t[2] = [float(x) for x in list(final_list_t[2])] * u.mas
        final_list_t[3] = [float(x) for x in list(final_list_t[3])] * u.mas
        final_list_t[4] = [float(x) for x in list(final_list_t[4])] * u.mas
        final_list_t[5] = [float(x) for x in list(final_list_t[5])] * u.mas
        final_list_t[7] = [float(x) for x in list(final_list_t[7])] * u.deg
        final_list_t[8] = [float(x) for x in list(final_list_t[8])] * u.deg
        final_list_t[9] = [float(x) for x in list(final_list_t[9])] * u.mas
        final_list_t[10] = [float(x) for x in list(final_list_t[10])] * u.mas
        out = QTable(final_list_t, names=('hip','hp_ref','d_ra*','d_dec','d_ra*_error','d_dec_error','f_d','position_angle',
                  'position_angle_error','separation_min','d_var','notes_flag','corr'),
            meta={'hip':'Hipparcos Catalogue identifier','hp_ref':'Adopted reference magnitudes for the object (mag)',
                  'd_ra*':'Element of the variability-induced motion in RA (mas)','d_dec':'Element of the variability-induced motion in Dec (mas)',
                  'd_ra*_error':'Standard error of the VIM element in RA (mas)','d_dec_error':'Standard error of the VIM element in Dec (mas)',
                  'f_d':'Significance of the VIM elements',
                  'position_angle':'position angle of the constant component of the binary with respect to the variable component, in degrees',
                  'position_angle_error':'Standard error of the position angle, in degrees',
                  'separation_min':'Lower limit for the separation of the binary (mas)',
                  'd_var':'Displacement of photocentre between the minimum and maximum luminosity of the system (mas)',
                  'notes_flag':'Flag indicating a note at the end of relevent volume','corr':'21 correlation coefficients in a non-linear coding'})

    elif flag == 'X':
        print(f'HIP {HIP} is in a stochastic solution.')
        data_list = [x.split(':',1)[1].lstrip() for x in text_list[1:5]]
        value_list = [re.compile("\s{3,}").split(x) for x in data_list]
        void_flag = np.array([len(x) == 1 for x in value_list])
        _ = [value_list[i].insert(0, np.nan) for i in np.where(void_flag)[0]]
        value_list_t = list(zip(*value_list))
        final_list = [value_list_t[0]]
        final_list_t = list(zip(*final_list))
        final_list_t[1] = [float(x) for x in list(final_list_t[1])] * u.mas
        final_list_t[2] = [float(x) for x in list(final_list_t[2])] * u.mas
        out = QTable(final_list_t, names=('hip','cosmic_error','cosmic_error_error','notes_flag'),
                    meta={'hip':'Hipparcos Catalogue identifier','cosmic_error':'Cosmic error, epsilon (mas)',
                          'cosmic_error_error':'Estimated standard error of the cosmic error, sigma_epsilon (mas)',
                          'notes_flag':'Flag indicating a note at the end of relevent volume'})
    
    else:
        print(f'HIP {HIP} is in a single star solution.')
        
    return out 

