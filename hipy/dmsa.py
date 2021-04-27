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
    # Convert the names of stars to HIP numbers, and get HIP numbers.
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
    
    # Get the data of Double and Multiple Systems Annex(dmsa) in the Hipparcos catalogue.
    url2 = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?dmId={HIP}'
    webpage = str(urllib.request.urlopen(url2).read())
    soup = bs4.BeautifulSoup(webpage, features='html.parser')
    text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\n'") 
    text_list = text.split('\\n')
    
    # DMSA/C: Component solutions, including 'optical' double stars and long-period binary and multiple systems.
    # A little complicated.
    if flag == 'C':
        print(f'HIP {HIP} is in a component solution.')
        # Remove notes('DN','GN','PN') from text_list.
        new_text_list = []
        notes_list = ['DN','GN','PN']
        for text in text_list:
            if not any(note in text for note in notes_list):
                new_text_list.append(text)
        # Separate the situations of 2,3,4 components by the length of new_text_list.
        # 2 components in one solution
        if len(new_text_list) == 92 or len(new_text_list) == 183:
            if len(new_text_list) == 183:
                del new_text_list[92:]
            data_list = [x.split(':',1)[1].lstrip() for x in new_text_list[1:38]] + [x.split(':',1)[1].lstrip() for x in new_text_list[49:75]]
            value_list = [re.compile("\s{3,}").split(x) for x in data_list]
            void_flag = np.array([len(x) == 1 for x in value_list])
            _ = [value_list[i].insert(0, np.nan) for i in np.where(void_flag)[0]]
            value_list_t = list(zip(*value_list))
            final_list = [value_list_t[0]]
            corr_list = new_text_list[88] + new_text_list[89] + new_text_list[90] + new_text_list[91]
            final_list = np.append(final_list,[[corr_list]],1)
            final_list_t = list(zip(*final_list))
            # Assign units to columns.
            index_mag = [14,15,16,17,18,19, 40,41,42,43,44,45]
            index_deg = [20,21,31, 46,47,57]
            index_mas = [22,25,26,27, 48,51,52,53]
            index_mas_yr = [23,24,28,29, 49,50,54,55]
            index_arcsec = [32, 58]
            index_deg_yr = [33, 59]
            index_arcsec_yr = [34, 60]
            for i in index_mag:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mag
            for i in index_deg:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.deg
            for i in index_mas:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas
            for i in index_mas_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas/u.yr
            for i in index_arcsec:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.arcsec
            for i in index_deg_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.deg/u.yr
            for i in index_arcsec_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.arcsec/u.yr

            out = QTable(final_list_t, 
                names=('ccdm','solution_identifier','solution_type','solution_source','solution_quality','notes_flag',
                       'n_s','n_c','n_p','n_r','comp','sequential_comp',
                       'comp_identifier','hip','hp_mag','hp_mag_error','bt_mag','bt_mag_error','vt_mag','vt_mag_error','ra','dec','parallax','pmra','pmdec',
                       'ra_error','dec_error','parallax_error','pmra_error','pmdec_error','refer_comp','position_angle','angular_separation',
                       'position_angle_change_rate','angular_separation_change_rate','sequential_record','status_flags',               
                       'sequential_comp_2',
                       'comp_identifier_2','hip_2','hp_mag_2','hp_mag_error_2','bt_mag_2','bt_mag_error_2','vt_mag_2','vt_mag_error_2','ra_2','dec_2','parallax_2','pmra_2','pmdec_2',
                       'ra_error_2','dec_error_2','parallax_error_2','pmra_error_2','pmdec_error_2','refer_comp_2','position_angle_2','angular_separation_2',
                       'position_angle_change_rate_2','angular_separation_change_rate_2','sequential_record_2','status_flags_2',               
                       'corr'),
                meta={'ccdm':'CCDM number','solution_identifier':'Solution identifier in this CCDM number (1 to n_s)',
                      'solution_type':'Type of solution (F,I,L)','solution_source':'Source of solution (C,F,N)',
                      'solution_quality':'Quality of solution (A,B,C,D)','notes_flag':'Flag indicating a note',
                      'n_s':'Number of solutions pertaining to the system','n_c':'Number of components in this solution',
                      'n_p':'Number of free parameters in this solution','n_r':'Number of correlation records in this solution',
                      'comp':"This field contains the word 'COMP'",'sequential_comp':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier':'Component identifier (A,B,C,...)','hip':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag':'Magnitude of component, Hp(mag)','hp_mag_error':'Standard error of Hp magnitude (mag)',
                      'bt_mag':'Magnitude of component, B_T(mag)','bt_mag_error':'Standard error of B_T magnitude (mag)',
                      'vt_mag':'Magnitude of component, V_T(mag)','vt_mag_error':'Standard error of V_T magnitude (mag)',
                      'ra':'right ascension (deg, J1991.25)','dec':'declination (deg, J1991.25)','parallax':'Trigonometric parallax (mas)',
                      'pmra':'proper motion in RA direction (mas/yr,J1991.25)','pmdec':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error':'Standard error of RA (mas, J1991.25)','dec_error':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error':'Standard error of parallax (mas)',
                      'pmra_error':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp':'Reference component for the relative data of the current component',
                      'position_angle':'Position angle relative to reference component (deg)','angular_separation':'Separation from reference component (arcsec)',
                      'position_angle_change_rate':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record':'Sequential record number for the reference component in refer_comp',
                      'status_flags':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'sequential_comp_2':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier_2':'Component identifier (A,B,C,...)','hip_2':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag_2':'Magnitude of component, Hp(mag)','hp_mag_error_2':'Standard error of Hp magnitude (mag)',
                      'bt_mag_2':'Magnitude of component, B_T(mag)','bt_mag_error_2':'Standard error of B_T magnitude (mag)',
                      'vt_mag_2':'Magnitude of component, V_T(mag)','vt_mag_error_2':'Standard error of V_T magnitude (mag)',
                      'ra_2':'right ascension (deg, J1991.25)','dec_2':'declination (deg, J1991.25)','parallax_2':'Trigonometric parallax (mas)',
                      'pmra_2':'proper motion in RA direction (mas/yr,J1991.25)','pmdec_2':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error_2':'Standard error of RA (mas, J1991.25)','dec_error_2':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error_2':'Standard error of parallax (mas)',
                      'pmra_error_2':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error_2':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp_2':'Reference component for the relative data of the current component',
                      'position_angle_2':'Position angle relative to reference component (deg)','angular_separation_2':'Separation from reference component (arcsec)',
                      'position_angle_change_rate_2':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate_2':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record_2':'Sequential record number for the reference component in refer_comp',
                      'status_flags_2':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'corr':'66 correlation coefficients in a non-linear coding'})
        # 3 components in one solution
        elif len(new_text_list) == 163 or len(new_text_list) == 254:
            if len(new_text_list) == 254:
                del new_text_list[163:]    
            data_list = [x.split(':',1)[1].lstrip() for x in new_text_list[1:38]] + [x.split(':',1)[1].lstrip() for x in new_text_list[49:75]]\
                      + [x.split(':',1)[1].lstrip() for x in new_text_list[86:112]]
            value_list = [re.compile("\s{3,}").split(x) for x in data_list]
            void_flag = np.array([len(x) == 1 for x in value_list])
            _ = [value_list[i].insert(0, np.nan) for i in np.where(void_flag)[0]]
            value_list_t = list(zip(*value_list))
            final_list = [value_list_t[0]]
            corr_list = new_text_list[125] + new_text_list[126] + new_text_list[127] + new_text_list[128]\
                      + new_text_list[142] + new_text_list[143] + new_text_list[144] + new_text_list[145]\
                      + new_text_list[159] + new_text_list[160] + new_text_list[161] + new_text_list[162]
            final_list = np.append(final_list,[[corr_list]],1)
            final_list_t = list(zip(*final_list))
            # Assign units to columns.
            index_mag = [14,15,16,17,18,19, 40,41,42,43,44,45, 66,67,68,69,70,71]
            index_deg = [20,21,31, 46,47,57, 72,73,83]
            index_mas = [22,25,26,27, 48,51,52,53, 74,77,78,79]
            index_mas_yr = [23,24,28,29, 49,50,54,55, 75,76,80,81]
            index_arcsec = [32, 58, 84]
            index_deg_yr = [33, 59, 85]
            index_arcsec_yr = [34, 60, 86]
            for i in index_mag:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mag
            for i in index_deg:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.deg
            for i in index_mas:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas
            for i in index_mas_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas/u.yr
            for i in index_arcsec:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.arcsec
            for i in index_deg_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.deg/u.yr
            for i in index_arcsec_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.arcsec/u.yr

            out = QTable(final_list_t,
                names=('ccdm','solution_identifier','solution_type','solution_source','solution_quality','notes_flag',
                       'n_s','n_c','n_p','n_r','comp','sequential_comp',
                       'comp_identifier','hip','hp_mag','hp_mag_error','bt_mag','bt_mag_error','vt_mag','vt_mag_error','ra','dec','parallax','pmra','pmdec',
                       'ra_error','dec_error','parallax_error','pmra_error','pmdec_error','refer_comp','position_angle','angular_separation',
                       'position_angle_change_rate','angular_separation_change_rate','sequential_record','status_flags',               
                       'sequential_comp_2',
                       'comp_identifier_2','hip_2','hp_mag_2','hp_mag_error_2','bt_mag_2','bt_mag_error_2','vt_mag_2','vt_mag_error_2','ra_2','dec_2','parallax_2','pmra_2','pmdec_2',
                       'ra_error_2','dec_error_2','parallax_error_2','pmra_error_2','pmdec_error_2','refer_comp_2','position_angle_2','angular_separation_2',
                       'position_angle_change_rate_2','angular_separation_change_rate_2','sequential_record_2','status_flags_2',
                       'sequential_comp_3',
                       'comp_identifier_3','hip_3','hp_mag_3','hp_mag_error_3','bt_mag_3','bt_mag_error_3','vt_mag_3','vt_mag_error_3','ra_3','dec_3','parallax_3','pmra_3','pmdec_3',
                       'ra_error_3','dec_error_3','parallax_error_3','pmra_error_3','pmdec_error_3','refer_comp_3','position_angle_3','angular_separation_3',
                       'position_angle_change_rate_3','angular_separation_change_rate_3','sequential_record_3','status_flags_3',
                       'corr'),
                meta={'ccdm':'CCDM number','solution_identifier':'Solution identifier in this CCDM number (1 to n_s)',
                      'solution_type':'Type of solution (F,I,L)','solution_source':'Source of solution (C,F,N)',
                      'solution_quality':'Quality of solution (A,B,C,D)','notes_flag':'Flag indicating a note',
                      'n_s':'Number of solutions pertaining to the system','n_c':'Number of components in this solution',
                      'n_p':'Number of free parameters in this solution','n_r':'Number of correlation records in this solution',
                      'comp':"This field contains the word 'COMP'",'sequential_comp':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier':'Component identifier (A,B,C,...)','hip':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag':'Magnitude of component, Hp(mag)','hp_mag_error':'Standard error of Hp magnitude (mag)',
                      'bt_mag':'Magnitude of component, B_T(mag)','bt_mag_error':'Standard error of B_T magnitude (mag)',
                      'vt_mag':'Magnitude of component, V_T(mag)','vt_mag_error':'Standard error of V_T magnitude (mag)',
                      'ra':'right ascension (deg, J1991.25)','dec':'declination (deg, J1991.25)','parallax':'Trigonometric parallax (mas)',
                      'pmra':'proper motion in RA direction (mas/yr,J1991.25)','pmdec':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error':'Standard error of RA (mas, J1991.25)','dec_error':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error':'Standard error of parallax (mas)',
                      'pmra_error':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp':'Reference component for the relative data of the current component',
                      'position_angle':'Position angle relative to reference component (deg)','angular_separation':'Separation from reference component (arcsec)',
                      'position_angle_change_rate':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record':'Sequential record number for the reference component in refer_comp',
                      'status_flags':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'sequential_comp_2':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier_2':'Component identifier (A,B,C,...)','hip_2':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag_2':'Magnitude of component, Hp(mag)','hp_mag_error_2':'Standard error of Hp magnitude (mag)',
                      'bt_mag_2':'Magnitude of component, B_T(mag)','bt_mag_error_2':'Standard error of B_T magnitude (mag)',
                      'vt_mag_2':'Magnitude of component, V_T(mag)','vt_mag_error_2':'Standard error of V_T magnitude (mag)',
                      'ra_2':'right ascension (deg, J1991.25)','dec_2':'declination (deg, J1991.25)','parallax_2':'Trigonometric parallax (mas)',
                      'pmra_2':'proper motion in RA direction (mas/yr,J1991.25)','pmdec_2':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error_2':'Standard error of RA (mas, J1991.25)','dec_error_2':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error_2':'Standard error of parallax (mas)',
                      'pmra_error_2':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error_2':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp_2':'Reference component for the relative data of the current component',
                      'position_angle_2':'Position angle relative to reference component (deg)','angular_separation_2':'Separation from reference component (arcsec)',
                      'position_angle_change_rate_2':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate_2':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record_2':'Sequential record number for the reference component in refer_comp',
                      'status_flags_2':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'sequential_comp_3':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier_3':'Component identifier (A,B,C,...)','hip_3':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag_3':'Magnitude of component, Hp(mag)','hp_mag_error_3':'Standard error of Hp magnitude (mag)',
                      'bt_mag_3':'Magnitude of component, B_T(mag)','bt_mag_error_3':'Standard error of B_T magnitude (mag)',
                      'vt_mag_3':'Magnitude of component, V_T(mag)','vt_mag_error_3':'Standard error of V_T magnitude (mag)',
                      'ra_3':'right ascension (deg, J1991.25)','dec_3':'declination (deg, J1991.25)','parallax_3':'Trigonometric parallax (mas)',
                      'pmra_3':'proper motion in RA direction (mas/yr,J1991.25)','pmdec_3':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error_3':'Standard error of RA (mas, J1991.25)','dec_error_3':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error_3':'Standard error of parallax (mas)',
                      'pmra_error_3':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error_3':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp_3':'Reference component for the relative data of the current component',
                      'position_angle_3':'Position angle relative to reference component (deg)','angular_separation_3':'Separation from reference component (arcsec)',
                      'position_angle_change_rate_3':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate_3':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record_3':'Sequential record number for the reference component in refer_comp',
                      'status_flags_3':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'corr':'153 correlation coefficients in a non-linear coding'})
        # 4 components in one solution
        elif len(new_text_list) == 234 or len(new_text_list) == 325:
            if len(new_text_list) == 325:
                del new_text_list[234:]
            data_list = [x.split(':',1)[1].lstrip() for x in new_text_list[1:38]] + [x.split(':',1)[1].lstrip() for x in new_text_list[49:75]]\
                      + [x.split(':',1)[1].lstrip() for x in new_text_list[86:112]] + [x.split(':',1)[1].lstrip() for x in new_text_list[123:149]]
            value_list = [re.compile("\s{3,}").split(x) for x in data_list]
            void_flag = np.array([len(x) == 1 for x in value_list])
            _ = [value_list[i].insert(0, np.nan) for i in np.where(void_flag)[0]]
            value_list_t = list(zip(*value_list))
            final_list = [value_list_t[0]]
            corr_list = new_text_list[162] + new_text_list[163] + new_text_list[164] + new_text_list[165]\
                      + new_text_list[179] + new_text_list[180] + new_text_list[181] + new_text_list[182]\
                      + new_text_list[196] + new_text_list[197] + new_text_list[198] + new_text_list[199]\
                      + new_text_list[213] + new_text_list[214] + new_text_list[215] + new_text_list[216]\
                      + new_text_list[230] + new_text_list[231] + new_text_list[232] + new_text_list[233]
            final_list = np.append(final_list,[[corr_list]],1)
            final_list_t = list(zip(*final_list))
            # Assign units to columns.
            index_mag = [14,15,16,17,18,19, 40,41,42,43,44,45, 66,67,68,69,70,71, 92,93,94,95,96,97]
            index_deg = [20,21,31, 46,47,57, 72,73,83, 98,99,109]
            index_mas = [22,25,26,27, 48,51,52,53, 74,77,78,79, 100,103,104,105]
            index_mas_yr = [23,24,28,29, 49,50,54,55, 75,76,80,81, 101,102,106,107]
            index_arcsec = [32, 58, 84, 110]
            index_deg_yr = [33, 59, 85, 111]
            index_arcsec_yr = [34, 60, 86, 112]
            for i in index_mag:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mag
            for i in index_deg:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.deg
            for i in index_mas:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas
            for i in index_mas_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas/u.yr
            for i in index_arcsec:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.arcsec
            for i in index_deg_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.deg/u.yr
            for i in index_arcsec_yr:
                final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.arcsec/u.yr

            out = QTable(final_list_t,
                names=('ccdm','solution_identifier','solution_type','solution_source','solution_quality','notes_flag',
                       'n_s','n_c','n_p','n_r','comp','sequential_comp',
                       'comp_identifier','hip','hp_mag','hp_mag_error','bt_mag','bt_mag_error','vt_mag','vt_mag_error','ra','dec','parallax','pmra','pmdec',
                       'ra_error','dec_error','parallax_error','pmra_error','pmdec_error','refer_comp','position_angle','angular_separation',
                       'position_angle_change_rate','angular_separation_change_rate','sequential_record','status_flags',               
                       'sequential_comp_2',
                       'comp_identifier_2','hip_2','hp_mag_2','hp_mag_error_2','bt_mag_2','bt_mag_error_2','vt_mag_2','vt_mag_error_2','ra_2','dec_2','parallax_2','pmra_2','pmdec_2',
                       'ra_error_2','dec_error_2','parallax_error_2','pmra_error_2','pmdec_error_2','refer_comp_2','position_angle_2','angular_separation_2',
                       'position_angle_change_rate_2','angular_separation_change_rate_2','sequential_record_2','status_flags_2',
                       'sequential_comp_3',
                       'comp_identifier_3','hip_3','hp_mag_3','hp_mag_error_3','bt_mag_3','bt_mag_error_3','vt_mag_3','vt_mag_error_3','ra_3','dec_3','parallax_3','pmra_3','pmdec_3',
                       'ra_error_3','dec_error_3','parallax_error_3','pmra_error_3','pmdec_error_3','refer_comp_3','position_angle_3','angular_separation_3',
                       'position_angle_change_rate_3','angular_separation_change_rate_3','sequential_record_3','status_flags_3',
                       'sequential_comp_4',
                       'comp_identifier_4','hip_4','hp_mag_4','hp_mag_error_4','bt_mag_4','bt_mag_error_4','vt_mag_4','vt_mag_error_4','ra_4','dec_4','parallax_4','pmra_4','pmdec_4',
                       'ra_error_4','dec_error_4','parallax_error_4','pmra_error_4','pmdec_error_4','refer_comp_4','position_angle_4','angular_separation_4',
                       'position_angle_change_rate_4','angular_separation_change_rate_4','sequential_record_4','status_flags_4',               
                       'corr'),
                meta={'ccdm':'CCDM number','solution_identifier':'Solution identifier in this CCDM number (1 to n_s)',
                      'solution_type':'Type of solution (F,I,L)','solution_source':'Source of solution (C,F,N)',
                      'solution_quality':'Quality of solution (A,B,C,D)','notes_flag':'Flag indicating a note',
                      'n_s':'Number of solutions pertaining to the system','n_c':'Number of components in this solution',
                      'n_p':'Number of free parameters in this solution','n_r':'Number of correlation records in this solution',
                      'comp':"This field contains the word 'COMP'",'sequential_comp':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier':'Component identifier (A,B,C,...)','hip':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag':'Magnitude of component, Hp(mag)','hp_mag_error':'Standard error of Hp magnitude (mag)',
                      'bt_mag':'Magnitude of component, B_T(mag)','bt_mag_error':'Standard error of B_T magnitude (mag)',
                      'vt_mag':'Magnitude of component, V_T(mag)','vt_mag_error':'Standard error of V_T magnitude (mag)',
                      'ra':'right ascension (deg, J1991.25)','dec':'declination (deg, J1991.25)','parallax':'Trigonometric parallax (mas)',
                      'pmra':'proper motion in RA direction (mas/yr,J1991.25)','pmdec':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error':'Standard error of RA (mas, J1991.25)','dec_error':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error':'Standard error of parallax (mas)',
                      'pmra_error':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp':'Reference component for the relative data of the current component',
                      'position_angle':'Position angle relative to reference component (deg)','angular_separation':'Separation from reference component (arcsec)',
                      'position_angle_change_rate':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record':'Sequential record number for the reference component in refer_comp',
                      'status_flags':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'sequential_comp_2':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier_2':'Component identifier (A,B,C,...)','hip_2':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag_2':'Magnitude of component, Hp(mag)','hp_mag_error_2':'Standard error of Hp magnitude (mag)',
                      'bt_mag_2':'Magnitude of component, B_T(mag)','bt_mag_error_2':'Standard error of B_T magnitude (mag)',
                      'vt_mag_2':'Magnitude of component, V_T(mag)','vt_mag_error_2':'Standard error of V_T magnitude (mag)',
                      'ra_2':'right ascension (deg, J1991.25)','dec_2':'declination (deg, J1991.25)','parallax_2':'Trigonometric parallax (mas)',
                      'pmra_2':'proper motion in RA direction (mas/yr,J1991.25)','pmdec_2':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error_2':'Standard error of RA (mas, J1991.25)','dec_error_2':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error_2':'Standard error of parallax (mas)',
                      'pmra_error_2':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error_2':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp_2':'Reference component for the relative data of the current component',
                      'position_angle_2':'Position angle relative to reference component (deg)','angular_separation_2':'Separation from reference component (arcsec)',
                      'position_angle_change_rate_2':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate_2':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record_2':'Sequential record number for the reference component in refer_comp',
                      'status_flags_2':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'sequential_comp_3':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier_3':'Component identifier (A,B,C,...)','hip_3':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag_3':'Magnitude of component, Hp(mag)','hp_mag_error_3':'Standard error of Hp magnitude (mag)',
                      'bt_mag_3':'Magnitude of component, B_T(mag)','bt_mag_error_3':'Standard error of B_T magnitude (mag)',
                      'vt_mag_3':'Magnitude of component, V_T(mag)','vt_mag_error_3':'Standard error of V_T magnitude (mag)',
                      'ra_3':'right ascension (deg, J1991.25)','dec_3':'declination (deg, J1991.25)','parallax_3':'Trigonometric parallax (mas)',
                      'pmra_3':'proper motion in RA direction (mas/yr,J1991.25)','pmdec_3':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error_3':'Standard error of RA (mas, J1991.25)','dec_error_3':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error_3':'Standard error of parallax (mas)',
                      'pmra_error_3':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error_3':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp_3':'Reference component for the relative data of the current component',
                      'position_angle_3':'Position angle relative to reference component (deg)','angular_separation_3':'Separation from reference component (arcsec)',
                      'position_angle_change_rate_3':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate_3':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record_3':'Sequential record number for the reference component in refer_comp',
                      'status_flags_3':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'sequential_comp_4':'Sequential component number in this solution (1 to n_c)',
                      'comp_identifier_4':'Component identifier (A,B,C,...)','hip_4':'Hipparcos Catalogue(HIP) identifier',
                      'hp_mag_4':'Magnitude of component, Hp(mag)','hp_mag_error_4':'Standard error of Hp magnitude (mag)',
                      'bt_mag_4':'Magnitude of component, B_T(mag)','bt_mag_error_4':'Standard error of B_T magnitude (mag)',
                      'vt_mag_4':'Magnitude of component, V_T(mag)','vt_mag_error_4':'Standard error of V_T magnitude (mag)',
                      'ra_4':'right ascension (deg, J1991.25)','dec_4':'declination (deg, J1991.25)','parallax_4':'Trigonometric parallax (mas)',
                      'pmra_4':'proper motion in RA direction (mas/yr,J1991.25)','pmdec_4':'proper motion in Dec direction (mas/yr,J1991.25)',
                      'ra_error_4':'Standard error of RA (mas, J1991.25)','dec_error_4':'Standard error of Dec (mas, J1991.25)',
                      'parallax_error_4':'Standard error of parallax (mas)',
                      'pmra_error_4':'Standard error of proper motion in RA direction (mas/yr)','pmdec_error_4':'Standard error of proper motion in Dec direction (mas/yr)',
                      'refer_comp_4':'Reference component for the relative data of the current component',
                      'position_angle_4':'Position angle relative to reference component (deg)','angular_separation_4':'Separation from reference component (arcsec)',
                      'position_angle_change_rate_4':'Rate of change of the position angle (deg/yr)',
                      'angular_separation_change_rate_4':'Rate of change of the angular separation (arcsec/yr)',
                      'sequential_record_4':'Sequential record number for the reference component in refer_comp',
                      'status_flags_4':'Status flags for Hp,ra,dec,parallax,pmra,pmdec: 1=estimated, 0=constrained to the value for the first component',

                      'corr':'276 correlation coefficients in a non-linear coding'})
        
        print(f'For more detailed information, please refer to https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?dmId={HIP}')
    
    # DMSA/G: Acceleration solutions. It lists apparently single(unresolved) stars, for which the motion appears to be significant non-linear.
    # They are probably 'astrometric binaries': either too close to be resolved (angular separation <= 0.1 arcsec), 
    # or with a companion too faint to be seen by Hipparcos.
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
        index_mas_yr2 = [1,2,3,4]
        index_mas_yr3 = [6,7,8,9]
        for i in index_mas_yr2:
            final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas/(u.yr)**2
        for i in index_mas_yr3:
            final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas/(u.yr)**3
        out = QTable(final_list_t, 
            names=('hip','g_ra*','g_dec','g_ra*_error','g_dec_error','f_g','gdot_ra*','gdot_dec',
                   'gdot_ra*_error','gdot_dec_error','f_gdot','notes_flag','n','corr'),
            meta={'hip':'Hipparcos Catalogue identifier','g_ra*':'Components in RA of the apparent acceleration of the potocentre at epoch J1991.25 (mas/yr^2)',
                 'g_dec':'Components in Dec of the apparent acceleration of the potocentre at epoch J1991.25 (mas/yr^2)',
                 'g_ra*_error':'Standard error of g_ra* (mas/yr^2)','g_dec_error':'Standard error of g_dec (mas/yr^2)','f_g':'Significance of the g terms',
                 'gdot_ra*':'Components in RA of the rate of change of the apparent acceleration of the photocentre (mas/yr^3)',
                 'gdot_dec':'Components in RA of the rate of change of the apparent acceleration of the photocentre (mas/yr^3)',
                 'gdot_ra*_error':'Standard error of g_dot_ra* (mas/yr^3)','gdot_dec_error':'Standard error of g_dot_dec (mas/yr^3)',
                  'f_gdot':'Significance of the gdot terms','notes_flag':'Flag indicating a note at the end of relevent volume',
                 'n':'Number of astrometric parameters, n = 7 or 9','corr':'(Up to) 36 correlation coefficients in a non-linear coding'})
        print(f'For more detailed information, please refer to https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?dmId={HIP}')

    # DMSA/O: Orbital solutions. It contains results for orbital binaries where the Hipparcos observations could be used to determine some or
    # all of the elements of the absolute Keplerian orbit of the system's photocentre.
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
        final_list_t[8] = [float(x) for x in list(final_list_t[8])] * u.d
        final_list_t[9] = [float(x) for x in list(final_list_t[9])] * u.d
        final_list_t[10] = [float(x) for x in list(final_list_t[10])] * u.mas
        final_list_t[11] = [float(x) for x in list(final_list_t[11])]
        index_deg = [5,6,7,12,13,14]
        for i in index_deg:
            final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.deg

        out = QTable(final_list_t, 
            names=('hip','p','t','a0','e','periastron_argument','i','node_position_angle',
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
        print(f'For more detailed information, please refer to https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?dmId={HIP}')

        
    # DMSA/V: VIM solutions, 'Variability-Induced' Movers. They are unresolved binaries in which one of the components is variable.
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
        final_list_t[7] = [float(x) for x in list(final_list_t[7])] * u.deg
        final_list_t[8] = [float(x) for x in list(final_list_t[8])] * u.deg
        index_mas = [2,3,4,5,9,10]
        for i in index_mas:
            final_list_t[i] = [float(x) for x in list(final_list_t[i])] * u.mas
        
        out = QTable(final_list_t, 
            names=('hip','hp_ref','d_ra*','d_dec','d_ra*_error','d_dec_error','f_d','position_angle',
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
        print(f'For more detailed information, please refer to https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?dmId={HIP}')

        
    # DMSA/X: Stochastic solutions. It finally lists the objects for which no solutions of the previous types could be found in reasonable
    # agreement with the standard errors of the Hipparcos observations. These objects could be double or multiple stars of unknown chracteristics,
    # or short-period astrometric binaries (< 3 years).
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
        out = QTable(final_list_t, 
                    names=('hip','cosmic_error','cosmic_error_error','notes_flag'),
                    meta={'hip':'Hipparcos Catalogue identifier','cosmic_error':'Cosmic error, epsilon (mas)',
                          'cosmic_error_error':'Estimated standard error of the cosmic error, sigma_epsilon (mas)',
                          'notes_flag':'Flag indicating a note at the end of relevent volume'})
        print(f'For more detailed information, please refer to https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?dmId={HIP}')
    

    else:
        print(f'HIP {HIP} is in a single star solution.')
        
    return out   