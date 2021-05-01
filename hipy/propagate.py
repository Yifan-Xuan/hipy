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
# Convert RA, Dec to Cartesian coordinates.
def radec2xyz(dec,ra):
    x = math.cos(dec)*math.cos(ra)
    y = math.cos(dec)*math.sin(ra)
    z = math.sin(dec)
    return [x,y,z]
# Convert Cartesian coordinates to RA, Dec.
def xyz2radec(x,y,z):
    dec = math.atan(z/math.sqrt(x*x+y*y))  # rad
    ra = math.atan(y/x)                    # rad
    if x < 0:
        ra = ra + math.pi                  # rad
    if (x >= 0) and (y < 0):
        ra = ra + 2*math.pi                # rad
    return [dec,ra]

d2r = math.pi/180
kpcmyr2auyr = 1e3*206265/1e6        # kpc/(10^6 yr) → au/yr
pc2au = 206265                      # pc → au
auyr2kms = 1.5*1e8/(365.25*24*3600) # au/yr → km/s
kpcmyr2kms = kpcmyr2auyr*auyr2kms   # kpc/(10^6 yr) → km/s
rv0 = 0  # km/s

# observational linear propagation
def propagate(star_name):
    
    # Convert the names of stars to HIP numbers, and get HIP numbers.
    if star_name.startswith('HIP'):
        HIP = int(star_name.split('HIP')[1])
    else:
        result_table = Simbad.query_objectids(star_name)
        line = list(filter(lambda x: 'HIP' in str(x), result_table))
        HIP = int(line[0][0].split('HIP')[1])

    print(f'### Propagation of HIP {HIP}')
    try:
        # Get RA, Dec, parallax, proper motion in RA and Dec at epoch J1991.25. 
        url = f'https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?noLinks=1&tabular=1&hipiId={HIP}'
        webpage = str(urllib.request.urlopen(url).read())
        soup = bs4.BeautifulSoup(webpage,'html.parser')
        text = soup.find(name='pre').get_text().lstrip("\\n").rstrip("\\r\\n\\r\\n\\r\\n'")
        text = text.split('\\r\\n\\r\\n')
        value_list = []
        value_list.append(text[0].split('\\n',2)[1].split('|'))
        ra0 = float(value_list[0][2])     # deg
        dec0 = float(value_list[0][3])    # deg
        plx0 = float(value_list[0][4])    # mas
        pmra0 = float(value_list[0][5])   # mas/yr
        pmdec0 = float(value_list[0][6])  # mas/yr
        ra0 = ra0*d2r                     # rad
        dec0 = dec0*d2r                   # rad
        # Get orbit number, source of abscissa, mid-epoch.
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
        orbit_number = data_list_t[0]
        source_absc = data_list_t[1]
        mid_epoch = [float(x) for x in data_list_t[10]]  # yr
        # Propagate observables to states.
        d0 = 1/plx0              # kpc
        r = radec2xyz(dec0,ra0) 
        x0 = r[0]*d0*1e3         # pc
        y0 = r[1]*d0*1e3         # pc
        z0 = r[2]*d0*1e3         # pc
        vdec = pmdec0*d0          # au/yr  # mas/yr*kpc = au/yr, the definition of parsec: 1 pc = 1 au / 1 arcsecond
        vra = pmra0*d0           # au/yr
        vp = math.sqrt(vra*vra + vdec*vdec)  # au/yr
        vr = rv0/auyr2kms    # au/yr
        vx = vr*math.cos(dec0)*math.cos(ra0) - vdec*math.sin(dec0)*math.cos(ra0) - vra*math.sin(ra0)  # au/yr          
        # note: vr is positive if the star is moving away from the Sun
        vy = vr*math.cos(dec0)*math.sin(ra0) - vdec*math.sin(dec0)*math.sin(ra0) + vra*math.cos(ra0)  # au/yr
        vz = vr*math.sin(dec0) + vdec*math.cos(dec0)  # au/yr
        x_list = []
        y_list = []
        z_list = []
        d_list = []
        ra_list = []
        dec_list = []
        plx_list = []
        for t in mid_epoch:
            x = x0 + vx*t/pc2au
            y = y0 + vy*t/pc2au
            z = z0 + vz*t/pc2au
            # Convert time-varying states to observables.
            ra = xyz2radec(x,y,z)[1]/d2r    # deg
            dec = xyz2radec(x,y,z)[0]/d2r   # deg
            d = math.sqrt(x*x+y*y+z*z)*1e-3 # kpc
            plx = 1.0/d                     # mas
            x_list.append(x)
            y_list.append(y)
            z_list.append(z)
            ra_list.append(ra)
            dec_list.append(dec)
            d_list.append(d)
            plx_list.append(plx)
        # velocity → proper motion
        # vequ = []
        pmra_list = []
        pmdec_list = []
        rv_list = []
        for j in range(len(mid_epoch)):
            rotz = np.matrix([[math.cos(ra_list[j]),math.sin(ra_list[j]),0.0],[-math.sin(ra_list[j]),math.cos(ra_list[j]),0.0],[0.0,0.0,1.0]])  # O-xyz → O-x'y'z'
            roty = np.matrix([[math.cos(dec_list[j]),0.0,math.sin(dec_list[j])],[0.0,1.0,0.0],[-math.sin(dec_list[j]),0.0,math.cos(dec_list[j])]])  # O-x'y'z' → O-x''y''z''
            c = np.array([vx,vy,vz]).reshape((-1,1))
            v = np.dot(np.dot(roty,rotz),c)
        #   vequ.append(v)
            pmra = v[1]/d_list[j]   # mas/yr
            pmdec = v[2]/d_list[j]  # mas/yr
            rv = v[0]*auyr2kms      # km/s
            pmra_list.append(pmra)
            pmdec_list.append(pmdec)
            rv_list.append(rv)
        final_list = [orbit_number,source_absc,ra_list,dec_list,plx_list,pmra_list,pmdec_list,rv_list,mid_epoch]
        final_list[2] = [float(x) for x in list(final_list[2])] * u.deg
        final_list[3] = [float(x) for x in list(final_list[3])] * u.deg
        final_list[4] = [float(x) for x in list(final_list[4])] * u.mas
        final_list[5] = [float(x) for x in list(final_list[5])] * u.mas/u.yr
        final_list[6] = [float(x) for x in list(final_list[6])] * u.mas/u.yr
        final_list[7] = [float(x) for x in list(final_list[7])] * u.km/u.s
        final_list[8] = [float(x) for x in list(final_list[8])] * u.yr
        out = QTable(final_list,
                    names=('orbit_number','source_absc','ra','dec','parallax','pmra','pmdec','rv','mid_epoch'),
                    meta={'orbit_number':'orbit number','source_absc':'source of abscissa (F or f if FAST data, N or n if NDAC data)',
                         'ra':'Right ascension in model','dec':'Declination in model','parallax':'parallax','pmra':'proper motion in RA','pmdec':'proper motion in Dec',
                         'rv':'Radial velocity','mid_epoch':'FAST/NDAC reference great-circle mid-epoch, in years relative to J1991.25(TT)'})
    except IndexError:
        print(f'The intermediate data of HIP {HIP} cannot be found.')
    return out