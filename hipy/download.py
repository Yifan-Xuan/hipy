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


# Construct functions to convert years(float) to Julian Date(float).
def cal2jd(yr,mn,dy):
    y = yr - 1
    m = mn + 12
#   date1 = 4.5+31*(10+12*1582);   # Last day of Julian calendar (1582.10.04 Noon)
#   date2 = 15.5+31*(10+12*1582);  # First day of Gregorian calendar (1582.10.15 Noon)
#   date = dy+31*(mn+12*yr)
    b = y//400 - y//100
#   b = math.floor(y/400) - math.floor(y/100)
    jd = math.floor(365.25*y) + math.floor(30.6001*(m+1)) + b + 1720996.5 + dy
    return jd
def yr2jd(yr):  
    iyr = math.floor(yr)
    jd0 = cal2jd(iyr,1,1)
    days = cal2jd(iyr+1,1,1)-jd0
    doy = (yr-iyr)*days+1
    jd1 = cal2jd(iyr,1,0) + doy
    return jd1

# This is from L. Lindegren's code to give the barycentric position of the earth in a given year between 1988 an 1993. 
def cearth(year):
    omega = 0.0172021240/86400
    e = 0.016714
    g0 = -0.04128
    au = 1.496e+11
    c = [[.540817e+09, -.334118e+11, -.145868e+12],
         [-.202315e+10, .133781e+12, -.306652e+11],
         [-.883589e+09, .580048e+11, -.132951e+11]]
    ret = []
    t  = (year - 1988.0) * 365.25 * 86400.0
    arg = omega * t + g0
    arg = arg + e * math.sin(arg)
    co = math.cos(arg)
    si = math.sin(arg)
    for i in range(3):
        ret.append((c[i][0] + c[i][1]*co + c[i][2]*si)/au)
    return ret

# Get five astrometric parameters of an HIP entry. 
def query(HIP):
    url = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?noLinks=1&tabular=1&hipiId={HIP}'
    webpage = str(urllib.request.urlopen(url).read())
    soup = bs4.BeautifulSoup(webpage,'html.parser')
    text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\r\\n\\r\\n\\r\\n'")
    text = text.split('\\r\\n\\r\\n')[0].split('\\n',1)[1].split('\\n',1)[0].split('|')
    p = [float(x) for x in text[2:7]]
    return p

def get_data(star_name,type):
    
    # Convert the names of stars to HIP numbers, and get HIP numbers.
    if star_name.startswith('HIP'):
        HIP = int(star_name.split('HIP')[1])
    else:
        result_table = Simbad.query_objectids(star_name)
        line = list(filter(lambda x: 'HIP' in str(x), result_table))
        HIP = int(line[0][0].split('HIP')[1])

    # Download 3 types of Hipparcos data: catalogue data (catalogue), intermediate data (intermediate), and epoch photometry data (epd).
    # Download Hipparcos catalogue data.
    if type == 'catalogue':
        try:
            print(f'### Query for catalogue_HIP {HIP}')
            url = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?hipId={HIP}'
            webpage = str(urllib.request.urlopen(url).read()) # Load the webpage
            # Parse webpage using https://www.crummy.com/software/BeautifulSoup/bs4/doc/#quick-start
            soup = bs4.BeautifulSoup(webpage, features='html.parser')
            # Get the text from the webpage. Note that the table is under the "pre" tag
            text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\n'") # remove extra '\\n' 
            text_list = text.split('\\n')
            # Now we have each row of the table. But still need to split each row

            # The first column and the first row are useless.
            # Notice that the last 4 lines are comments, thus useless.
            # `data_list` is just useful contents
            data_list = [x.split(':',1) for x in text_list[1:79]] # first split by ":"
            data_list_t = list(zip(*data_list)) # transpose the nested list
            data_list_t[1] = [x.lstrip(' ') for x in data_list_t[1]] # strip redundant whitespaces from left-hand side

            value_list = [re.compile("\s{3,}").split(x) for x in data_list_t[1][0:73]] + [re.compile("\s{6,}").split(x) for x in data_list_t[1][73:]]
            # split the string after ":" by requiring at least 3 or 6 whitespaces. 
            void_flag = np.array([len(x) == 1 for x in value_list]) # Some rows do not have value, only have field description. This array is a mask for them.
            _ = [value_list[i].insert(0, np.nan) for i in np.where(void_flag)[0]] # For rows without value, assign np.nan
            value_list_t = list(zip(*value_list)) # Transpose the value_list
            final_list = [value_list_t[0]]
            final_list_t = list(zip(*final_list))
            
            # Assign units to columns.
            index_mag = [5,32,33,34,35,37,38,40,41,44,45,46,49,50,66,67,75]
            index_deg = [8,9,63]
            index_mas = [11,14,15,16]
            index_mas_yr = [12,13,17,18]
            index_float = [19,20,21,22,23,24,25,26,27,28,29,30]
            for i in index_mag:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mag
            for i in index_deg:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.deg
            for i in index_mas:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas
            for i in index_mas_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas/u.yr
            for i in index_float:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])]
            final_list_t[51] = [float(x) for x in list(final_list_t[51])] * u.d
            final_list_t[64] = [float(x) for x in list(final_list_t[64])] * u.arcsec
            final_list_t[65] = [float(x) for x in list(final_list_t[65])] * u.arcsec
            
            # Assign names and descriptions to columns.
            out = QTable(final_list_t,
                 names=('catalogue','hip','proximity_flag','ra_hms','dec_dms','v_mag','coarse_varflag','v_mag_source',
                       'ra','dec','astrometry_flag','parallax','pmra','pmdec','ra_error','dec_error','parallax_error','pmra_error','pmdec_error',                        
                       'ra_dec_corr','ra_parallax_corr','dec_parallax_corr','ra_pmra_corr','dec_pmra_corr','parallax_pmra_corr','ra_pmdec_corr','dec_pmdec_corr','parallax_pmdec_corr','pmra_pmdec_corr',                      
                        'f1','f2','hip_number','bt_mag','bt_mag_error','vt_mag','vt_mag_error','bt_vt_flag','b_v','b_v_error','b_v_source','v_i','v_i_error','v_i_source','colour_indices_flag',             
                        'hp_mag','hp_mag_error','hp_scatter','hp_number','phot_flag','hp_max','hp_min','var_period','var_type','var_tables_flag','var_light_curves_flag',                        
                       'ccdm','ccdm_historical_status','n_catalogue_entries','n_components','dmsa_flag','astrometric_source_flag','solution_quality','component_identifiers',
                        'position_angle','angular_separation','angular_separation_error','delta_hp','delta_hp_error',
                        'survey_flag','chart_flag','notes_flag','hd','bd','cod','cpd','v_i_mag_red','sp_type','sp_type_source'),

                 meta={'catalogue':'Catalogue (H = Hipparcos, T = Tycho)', 'hip':'Hipparcos Catalogue(HIP) identifier', 'proximity_flag':'Proximity flag',
                       'ra_hms':'The approximate right ascension in conventional sexagesimal units with truncated precision and within the ICRS reference system, for epoch J1991.25', 
                       'dec_dms':'The approximate declination in conventional sexagesimal units with truncated precision and within the ICRS reference system, for epoch J1991.25',
                      'v_mag':'The magnitude, V, in the Johnson UBV photometric system','coarse_varflag':'Coarse variability flag', 'v_mag_source':'Source of magnitude identifier',
                      'ra':'right ascension in degrees (J1991.25, ICRS)','dec':'declination in degrees (J1991.25, ICRS)','astrometry_flag':'Reference flag for astrometric parameters of double and multiple systems',
                      'parallax':'Trigonometric parallax (mas)','pmra':'proper motion in RA direction (mas/yr,J1991.25,ICRS)','pmdec':'proper motion in Dec direction (mas/yr,J1991.25,ICRS)',
                       'ra_error':'Standard error of RA (mas, J1991.25)','dec_error':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error':'Standard error of parallax (mas)','pmra_error':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error':'Standard error of proper motion in Dec direction (mas/yr)',
                      'ra_dec_corr':'correlation coefficient between RA and Dec','ra_parallax_corr':'correlation coefficient between RA and parallax','dec_parallax_corr':'correlation coefficient between declination and parallax',          
                       'ra_pmra_corr':'correlation coefficient between RA and pmra','dec_pmra_corr':'correlation coefficient between Dec and pmra','parallax_pmra_corr':'correlation coefficient between parallax and pmra',
                       'ra_pmdec_corr':'correlation coefficient between RA and pmdec','dec_pmdec_corr':'correlation coefficient between Dec and pmdec','parallax_pmdec_corr':'correlation coefficient between parallax and pmdec','pmra_pmdec_corr':'correlation coefficient between pmra and pmdec',                      
                       'f1':'The percentage of rejected data','f2':'Goodness-of-fit statistic','hip_number':'Hipparcos Catalogue(HIP) identifier','bt_mag':'BT, Mean magnitude in the Tycho photometric system (mag)','bt_mag_error':'Standard error of the BT magnitude (mag)',
                       'vt_mag':'VT, Mean magnitude in the Tycho photometric system (mag)','vt_mag_error':'Standard error of the VT magnitude (mag)','bt_vt_flag':'Reference flag for BT and VT','b_v':'Colour index, B-V (mag)','b_v_error':'Standard error of the colour index (mag)','b_v_source':'Source of B-V',
                       'v_i':'Colour index, V-I (mag)','v_i_error':'Standard error of the colour index (mag)','v_i_source':'Source of the colour index, V-I','colour_indices_flag':'Reference flag for colour indices',            
                       'hp_mag':'Median magnitude in the Hipparcos photometric system, Hp (mag)','hp_mag_error':'Standard error of the median Hp magnitude (mag)','hp_scatter':'Scatter of the Hp observations, s (mag)','hp_number':'Number of Hp observations, N','phot_flag':'Reference flag for the photometry parameters',
                       'hp_max':'Mag at max, Hp (5th percentile)','hp_min':'Mag at min, Hp (95th percentile)','var_period':'Variability period from Hipparcos observations (days)','var_type':'Type of variability','var_tables_flag':'Variability annex flag: tabular data','var_light_curves_flag':'Variability annex flag: light curves',                       
                       'ccdm':'the identifier of the Catalogue of Component of Double and Multiple Stars','ccdm_historical_status':'Historical status of the CCDM identifier','n_catalogue_entries':'Number of separate catalogue entries with the same CCDM identifier','n_components':'Number of components into which the entry was resolved',
                       'dmsa_flag':'Double and multiple systems annex flag','astrometric_source_flag':'Source of the absolute astrometry','solution_quality':'Solution quality flag','component_identifiers':'Component designation for parameters below', 
                       'position_angle':'Position angle between the components, rounded (degrees,J1991.25)','angular_separation':'Angular separation between the components, rounded (arcsec,J1991.25)','angular_separation_error':'Standard error of the angular separation (arcsec)','delta_hp':'Magnitude difference of components (mag)','delta_hp_error':'Standard error of the magnitude difference (mag)',
                       'survey_flag':'Flag indicating survey star','chart_flag':'Flag indicating identification chart','notes_flag':'Flag indicating a note given at the end of the volumes','hd':'HD/HDE/HDEC identifier','bd':'DM identifier (BD)','cod':'DM identifier (CoD)','cpd':'DM identifier (CPD)',
                       'v_i_mag_red':'V-I (mag) used for photometric processing','sp_type':'Spectral type','sp_type_source':'Source of spectral type'})
            
            # Give the solution type of this HIP entry.
            if final_list_t[59][0] == 'C':
                print(f'HIP {HIP} is in a component solution.')
            elif final_list_t[59][0] == 'G':
                print(f'HIP {HIP} is in an acceleration solution.')
            elif final_list_t[59][0] == 'O':
                print(f'HIP {HIP} is in an orbital solution.')
            elif final_list_t[59][0] == 'V':
                print(f'HIP {HIP} is in a solution of variability-induced movers.')
            elif final_list_t[59][0] == 'X':
                print(f'HIP {HIP} is in a stochastic solution.')
            else:
                print(f'HIP {HIP} is in a single star solution.')
            
            print(f'For more detailed information, please refer to https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?hipId={HIP}')
        # A small minority of HIP catalogue data cannot be found.
        except IndexError:
            print(f'The catalogue of HIP {HIP} cannot be found.')

    # Download Hipparcos intermediate data.
    elif type == 'intermediate':
        try:
            print(f'### Query for intermediate_HIP {HIP}')
            url = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?noLinks=1&tabular=1&hipiId={HIP}'
            webpage = str(urllib.request.urlopen(url).read())
            soup = bs4.BeautifulSoup(webpage,'html.parser')
            text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\r\\n\\r\\n\\r\\n'")
            text = text.split('\\r\\n\\r\\n')
            text_list = []
            text_list.append(text[0].split('\\n',3)[3])
            for line in text[1:]:
                text_list.append(line)
            data_list = []
            for row in text_list:
                row = row.replace('\\n','|').split('|')
                if ('F' in row) or ('f' in row):
                    data = row[0:9]
                    data.append(np.nan) if row[9] == '     ' else data.append(row[9])
                    for x in row[11:14]:
                        data.append(x)
                elif ('N' in row) or ('n' in row):
                    data = row[0:9]
                    data.append(np.nan) if row[9] == '     ' else data.append(row[9])
                    for x in row[14:]:
                        data.append(x)
                data_list.append(data)
            data_list_t = list(map(list, zip(*data_list)))
            index_float = [2,3,7,8,9,10]
            for i in index_float:
                data_list_t[i] = [float(x) for x in list(data_list_t[i])]

            epoch_time = [x+1991.25 for x in list(data_list_t[10])]
            data_list_t = np.array(data_list_t).transpose()
            data_list_t = np.insert(data_list_t,11,epoch_time,axis=1)

            # Convert Julian years into Barycentric Julian Days.
            jd = []
            for yr in epoch_time:  
                jd.append(yr2jd(yr))
            data_list_t = np.insert(data_list_t,12,jd,axis=1)

            # Convert abscissa residuals into RA residuals and Dec residuals.
            ra_residual = [float(data_list_t[i][2])*float(data_list_t[i][7]) for i in range(np.size(data_list_t,0))]
            dec_residual = [float(data_list_t[i][3])*float(data_list_t[i][7]) for i in range(np.size(data_list_t,0))]
            data_list_t = np.insert(data_list_t,8,ra_residual,axis=1)
            data_list_t = np.insert(data_list_t,9,dec_residual,axis=1)

            # Convert standard errors of abscissa into standard errors of RA and Dec.
            ra_error = [float(data_list_t[i][2])*float(data_list_t[i][10]) for i in range(np.size(data_list_t,0))]
            dec_error = [float(data_list_t[i][3])*float(data_list_t[i][10]) for i in range(np.size(data_list_t,0))]
            data_list_t = np.insert(data_list_t,11,ra_error,axis=1)
            data_list_t = np.insert(data_list_t,12,dec_error,axis=1)

            # Get the correlation coefficients between RA and Dec.
            ra_dec_corr = [float(data_list_t[i][2])*float(data_list_t[i][3]) for i in range(np.size(data_list_t,0))]
            data_list_t = np.insert(data_list_t,14,ra_dec_corr,axis=1)

            # Get the observational values of ra and dec in corresponding FAST/NDAC great-circle epoch time.
            d2r = math.pi/180
            p = query(HIP)

            i = 0
            obs_dec_list = []
            obs_ra_list = []
            for year in epoch_time:
                ret = cearth(year)
                pa = ret[0]*math.sin(p[0]*d2r) - ret[1]*math.cos(p[0]*d2r)
                pd = (ret[0]*math.cos(p[0]*d2r) + ret[1]*math.sin(p[0]*d2r))*math.sin(p[1]*d2r) - ret[2]*math.cos(p[1]*d2r)
                model_dec = p[1] + (float(data_list_t[i][15])*p[4] + pd*p[2])/3600/1000
                obs_dec = model_dec + dec_residual[i]/3600/1000
                obs_ra = p[0] + (ra_residual[i] + (float(data_list_t[i][15])*p[3] + p[2]*pa)/math.cos(model_dec*d2r))/3600/1000 
                obs_dec_list.append(obs_dec)
                obs_ra_list.append(obs_ra)
                i = i + 1 
            data_list_t = np.insert(data_list_t,2,obs_ra_list,axis=1)
            data_list_t = np.insert(data_list_t,3,obs_dec_list,axis=1)

            data_list_tt = list(zip(*data_list_t))
            # Assign units to columns.
            index_deg = [2,3,20,21]
            index_float = [4,5,6,7,8,15,16]
            index_mas = [9,10,11,12,13,14]
            for i in index_deg:
                data_list_tt[i] = [float(x) for x in list(data_list_tt[i])] * u.deg
            for i in index_float:
                data_list_tt[i] = [float(x) for x in list(data_list_tt[i])]
            for i in index_mas:
                data_list_tt[i] = [float(x) for x in list(data_list_tt[i])] * u.mas
            data_list_tt[17] = [float(x) for x in list(data_list_tt[17])] * u.yr
            data_list_tt[18] = [float(x) for x in list(data_list_tt[18])] * u.yr
            data_list_tt[19] = [float(x) for x in list(data_list_tt[19])] * u.d

            # Assign names and descriptions to columns.
            out = QTable(data_list_tt,
                names=('orbit_number','source_absc','obs_ra','obs_dec','absc/ra','absc/dec','absc/parallax','absc/pmra','absc/pmdec','absc_residual','ra_residual',
                       'dec_residual','absc_error','ra_error','dec_error','absc_corr','ra_dec_corr','ref_great-circle_mid-epoch (yr)',
                        'ref_great-circle_epoch_time (yr)','ref_great-circle_epoch_time (bjd)','great-circle_pole_ra','dec_great-circle_pole_dec'),
                meta={'orbit_number':'orbit number','source_absc':'source of abscissa (F or f if FAST data, N or n if NDAC data)',
                      'obs_ra':'observational value of RA in this FAST/NDAC great-circle epoch time (deg)','obs_dec':'observational value of Dec in this FAST/NDAC great-circle epoch time (deg)',
                      'absc/ra':'abscissa partial derivative with respect to RA','absc/dec':'abscissa partial derivative with respect to Dec',
                      'absc/parallax':'abscissa partial derivative with respect to parallax','absc/pmra':'abscissa partial derivative with respect to proper motion in RA direction',
                      'absc/pmdec':'abscissa partial derivative with respect to proper motion in Dec direction',
                      'absc_residual':'abscissa residual','ra_residual':'residual of right ascension','dec_residual':'residual of declination',
                      'absc_error':'standard error of the abscissa','ra_error':'standard error of right ascension','dec_error':'standard error of declination',
                      'absc_corr':'correlation coefficient between FAST and NDAC abscissae','ra_dec_corr':'correlation coefficient between right ascension and declination',
                      'ref_great-circle_mid-epoch (yr)':'FAST/NDAC reference great-circle mid-epoch, in years relative to J1991.25(TT)',
                      'ref_great-circle_epoch_time (yr)':'the epoch time of the FAST/NDAC reference great-circle, in years',
                      'ref_great-circle_epoch_time (bjd)':'the epoch time of the FAST/NDAC reference great-circle, in Barycentric Julian Days',
                      'great-circle_pole_ra':'right ascension within ICRS of the FAST/NDAC reference great-circle pole',
                      'great-circle_pole_dec':'declination within ICRS of the FAST/NDAC reference great-circle pole'})            
            print(f'For more detailed information, please refer to https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?noLinks=1&tabular=1&hipiId={HIP}')
        # A small minority of HIP intermediate data cannot be found.
        except IndexError:
            print(f'The intermediate data of HIP {HIP} cannot be found.')

    # Download Hipparcos epoch photometry data.
    elif type == 'epd':
        try:
            print(f'### Query for epd_HIP {HIP}')
            url = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?hipepId={HIP}'
            webpage = str(urllib.request.urlopen(url).read())
            soup = bs4.BeautifulSoup(webpage,'html.parser')
            text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\r\\n'")
            text = text.split('\\n',17)[17]
            text_list = text.split('\\r\\n')
            data_list = [x.split('|') for x in text_list]
            data_list_t = list(map(list, zip(*data_list)))
            # Assign units to columns.
            data_list_t[0] = [float(x)+2440000 for x in list(data_list_t[0])] * u.d
            data_list_t[1] = [float(x) for x in list(data_list_t[1])] * u.mag
            data_list_t[2] = [float(x) for x in list(data_list_t[2])] * u.mag
            # Assign names and descriptions to columns.
            out = QTable(data_list_t,
                names=('obs_epoch','hp','hp_error','quality_flag'),
                meta={'obs_epoch':'observation epoch, in Barycentric Julian Days',
                      'hp':'calibrated Hp magnitude for this transit','hp_error':'estimated standard error of Hp magnitude','quality_flag':'quality flag, from bit 0 to bit 8'})
            
            print(f'For more detailed information, please refer to https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?hipepId={HIP}')
        # A small minority of HIP epoch photometry data cannot be found.
        except IndexError:
            print(f'The epoch photometry data of HIP {HIP} cannot be found.')
            
    else:
        print('No such type of data exits. Please check your input type.')

    return out
