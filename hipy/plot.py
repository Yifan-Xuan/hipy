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

def plot(star_name,consortium):

    # Convert the names of stars to HIP numbers, and get HIP numbers.
    if star_name.startswith('HIP'):
        HIP = int(star_name.split('HIP')[1])
    else:
        result_table = Simbad.query_objectids(star_name)
        line = list(filter(lambda x: 'HIP' in str(x), result_table))
        HIP = int(line[0][0].split('HIP')[1])

    try:
        # Get the barycentric and observed lines(solution) of motion for Hipparcos data.
        d2r = math.pi/180
        p = query(HIP)
        t0 = 1991.25
        year = np.linspace(1989.8,1993.3,341,endpoint=True)
        barycentric_ddec_list = []
        barycentric_dra_list = []
        model_ddec_list = []
        model_dra_list = []
        for t in year:
            ret = cearth(t)
            pa = ret[0]*math.sin(p[0]*d2r) - ret[1]*math.cos(p[0]*d2r)
            pd = (ret[0]*math.cos(p[0]*d2r) + ret[1]*math.sin(p[0]*d2r))*math.sin(p[1]*d2r) - ret[2]*math.cos(p[1]*d2r)
            barycentric_ddec = (t-t0)*p[4]
            barycentric_dra = (t-t0)*p[3]
            barycentric_ddec_list.append(barycentric_ddec)
            barycentric_dra_list.append(barycentric_dra)
            model_ddec = (t-t0)*p[4] + p[2]*pd
            model_dra = (t-t0)*p[3] + p[2]*pa
            model_ddec_list.append(model_ddec)
            model_dra_list.append(model_dra)

        # Get text from webpage. 
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
        data_list_t = np.array(list(map(list,zip(*data_list))))
        # Derive residuals and standard errors.
        source_absc = data_list_t[1].tolist()
        mid_epoch = [float(x) for x in data_list_t[10].tolist()]
        ra_residual = [float(data_list_t[2][i])*float(data_list_t[7][i]) for i in range(np.size(data_list_t,1))]
        dec_residual = [float(data_list_t[3][i])*float(data_list_t[7][i]) for i in range(np.size(data_list_t,1))]
        ra_error = [float(data_list_t[2][i])*float(data_list_t[8][i]) for i in range(np.size(data_list_t,1))]
        dec_error = [float(data_list_t[3][i])*float(data_list_t[8][i]) for i in range(np.size(data_list_t,1))]

        # Separate source of abscissa, mid_epoch, RA/Dec residual, RA/Dec standard error of the FAST consortium and the NDAC consortium.
        f_source_absc = []
        f_mid_epoch = []
        f_ra_residual = []
        f_dec_residual = []
        f_ra_error = []
        f_dec_error = []
        n_source_absc = []
        n_mid_epoch = []
        n_ra_residual = []
        n_dec_residual = []
        n_ra_error = []
        n_dec_error = []
        for i in range(len(source_absc)):
            if source_absc[i] == 'F' or source_absc[i] == 'f':
                f_source_absc.append(source_absc[i])
                f_mid_epoch.append(mid_epoch[i])
                f_ra_residual.append(ra_residual[i])
                f_dec_residual.append(dec_residual[i])
                f_ra_error.append(ra_error[i])
                f_dec_error.append(dec_error[i])
            elif source_absc[i] == 'N' or source_absc[i] == 'n':
                n_source_absc.append(source_absc[i])
                n_mid_epoch.append(mid_epoch[i])
                n_ra_residual.append(ra_residual[i])
                n_dec_residual.append(dec_residual[i])
                n_ra_error.append(ra_error[i])
                n_dec_error.append(dec_error[i])
        # Get the fitted and rejected points on the observed line(solution).
        f_fit_ddec_list = []
        f_fit_dra_list = []
        n_fit_ddec_list = []
        n_fit_dra_list = []
        for f in f_mid_epoch:
            ret = cearth(f+t0)
            pa = ret[0]*math.sin(p[0]*d2r) - ret[1]*math.cos(p[0]*d2r)
            pd = (ret[0]*math.cos(p[0]*d2r) + ret[1]*math.sin(p[0]*d2r))*math.sin(p[1]*d2r) - ret[2]*math.cos(p[1]*d2r)
            f_fit_ddec = f*p[4] + p[2]*pd
            f_fit_dra = f*p[3] + p[2]*pa
            f_fit_ddec_list.append(f_fit_ddec)
            f_fit_dra_list.append(f_fit_dra)
        for n in n_mid_epoch:
            ret = cearth(n+t0)
            pa = ret[0]*math.sin(p[0]*d2r) - ret[1]*math.cos(p[0]*d2r)
            pd = (ret[0]*math.cos(p[0]*d2r) + ret[1]*math.sin(p[0]*d2r))*math.sin(p[1]*d2r) - ret[2]*math.cos(p[1]*d2r)
            n_fit_ddec = n*p[4] + p[2]*pd
            n_fit_dra = n*p[3] + p[2]*pa
            n_fit_ddec_list.append(n_fit_ddec)
            n_fit_dra_list.append(n_fit_dra)
        # Plot the figure.
        if consortium == 'F' or consortium == 'f' or consortium == 'FAST':
            try:
                plt.figure(figsize=(10,10))
                for i in range(len(f_mid_epoch)):
                    plt.plot([f_fit_dra_list[i],f_fit_dra_list[i]+f_ra_residual[i]],[f_fit_ddec_list[i],f_fit_ddec_list[i]+f_dec_residual[i]],color='blue')
                    plt.plot([f_fit_dra_list[i]+f_ra_residual[i]+f_dec_error[i],f_fit_dra_list[i]+f_ra_residual[i]-f_dec_error[i]],
                             [f_fit_ddec_list[i]+f_dec_residual[i]-f_ra_error[i],f_fit_ddec_list[i]+f_dec_residual[i]+f_ra_error[i]],color='blue')
                    if f_source_absc[i] == 'F':
                        type1 = plt.scatter(f_fit_dra_list[i],f_fit_ddec_list[i],s=10,color='blue')
                    elif f_source_absc[i] == 'f':
                        type2 = plt.scatter(f_fit_dra_list[i],f_fit_ddec_list[i],s=100,facecolors='none',edgecolors='blue')
            except IndexError:
                print(f'The intermediate data of HIP {HIP} has no source of abscissa from F or f.')
            finally:
                plt.title(f'HIP {HIP}')
                plt.ylabel(r'$\Delta\delta$'+' [mas]')
                plt.xlabel(r'$\Delta\alpha$'+r'$ cos(\delta)$'+' [mas]')
                type5 = plt.plot(barycentric_dra_list,barycentric_ddec_list,color='violet',label='Barycentric motion')
                type6 = plt.plot(model_dra_list,model_ddec_list,color='lime',label='Solution')
                try:
                    plt.legend((type1,type2,type5[0],type6[0]),('FAST Data Fitted','FAST Data Rejected','Barycentric motion','Solution'))
                except NameError:
                    plt.legend((type1,type5[0],type6[0]),('FAST Data Fitted','Barycentric motion','Solution'))
                plt.show()
        elif consortium == 'N' or consortium == 'n' or consortium == 'NDAC':
            try:
                plt.figure(figsize=(10,10))
                for j in range(len(n_mid_epoch)):
                    plt.plot([n_fit_dra_list[j],n_fit_dra_list[j]+n_ra_residual[j]],[n_fit_ddec_list[j],n_fit_ddec_list[j]+n_dec_residual[j]],color='red')
                    plt.plot([n_fit_dra_list[j]+n_ra_residual[j]+n_dec_error[j],n_fit_dra_list[j]+n_ra_residual[j]-n_dec_error[j]],
                             [n_fit_ddec_list[j]+n_dec_residual[j]-n_ra_error[j],n_fit_ddec_list[j]+n_dec_residual[j]+n_ra_error[j]],color='red')
                    if n_source_absc[j] == 'N':
                        type3 = plt.scatter(n_fit_dra_list[j],n_fit_ddec_list[j],s=10,color='red')
                    elif n_source_absc[j] == 'n':
                        type4 = plt.scatter(n_fit_dra_list[j],n_fit_ddec_list[j],s=100,facecolors='none',edgecolors='red')
            except IndexError:
                print(f'The intermediate data of HIP {HIP} has no source of abscissa from N or n.')
            finally:
                plt.title(f'HIP {HIP}')
                plt.ylabel(r'$\Delta\delta$'+' [mas]')
                plt.xlabel(r'$\Delta\alpha$'+r'$ cos(\delta)$'+' [mas]')
                type5 = plt.plot(barycentric_dra_list,barycentric_ddec_list,color='violet',label='Barycentric motion')
                type6 = plt.plot(model_dra_list,model_ddec_list,color='lime',label='Solution')
                try:
                    plt.legend((type3,type4,type5[0],type6[0]),('NDAC Data Fitted','NDAC Data Rejected','Barycentric motion','Solution'))
                except NameError:
                    plt.legend((type3,type5[0],type6[0]),('NDAC Data Fitted','Barycentric motion','Solution'))
                plt.show()
        elif consortium == 'Both' or consortium == 'both':
            try:
                plt.figure(figsize=(10,10))
                for i in range(len(f_mid_epoch)):
                    plt.plot([f_fit_dra_list[i],f_fit_dra_list[i]+f_ra_residual[i]],[f_fit_ddec_list[i],f_fit_ddec_list[i]+f_dec_residual[i]],color='blue')
                    plt.plot([f_fit_dra_list[i]+f_ra_residual[i]+f_dec_error[i],f_fit_dra_list[i]+f_ra_residual[i]-f_dec_error[i]],
                             [f_fit_ddec_list[i]+f_dec_residual[i]-f_ra_error[i],f_fit_ddec_list[i]+f_dec_residual[i]+f_ra_error[i]],color='blue')
                    if f_source_absc[i] == 'F':
                        type1 = plt.scatter(f_fit_dra_list[i],f_fit_ddec_list[i],s=10,color='blue')
                    elif f_source_absc[i] == 'f':
                        type2 = plt.scatter(f_fit_dra_list[i],f_fit_ddec_list[i],s=100,facecolors='none',edgecolors='blue')
            except IndexError:
                print(f'The intermediate data of HIP {HIP} has no source of abscissa from F or f.')
            try:
                for j in range(len(n_mid_epoch)):
                    plt.plot([n_fit_dra_list[j],n_fit_dra_list[j]+n_ra_residual[j]],[n_fit_ddec_list[j],n_fit_ddec_list[j]+n_dec_residual[j]],color='red')
                    plt.plot([n_fit_dra_list[j]+n_ra_residual[j]+n_dec_error[j],n_fit_dra_list[j]+n_ra_residual[j]-n_dec_error[j]],
                             [n_fit_ddec_list[j]+n_dec_residual[j]-n_ra_error[j],n_fit_ddec_list[j]+n_dec_residual[j]+n_ra_error[j]],color='red')
                    if n_source_absc[j] == 'N':
                        type3 = plt.scatter(n_fit_dra_list[j],n_fit_ddec_list[j],s=10,color='red')
                    elif n_source_absc[j] == 'n':
                        type4 = plt.scatter(n_fit_dra_list[j],n_fit_ddec_list[j],s=100,facecolors='none',edgecolors='red')
            except IndexError:
                print(f'The intermediate data of HIP {HIP} has no source of abscissa from N or n.')    
            finally:
                plt.title(f'HIP {HIP}')
                plt.ylabel(r'$\Delta\delta$'+' [mas]')
                plt.xlabel(r'$\Delta\alpha$'+r'$ cos(\delta)$'+' [mas]')
                type5 = plt.plot(barycentric_dra_list,barycentric_ddec_list,color='violet',label='Barycentric motion')
                type6 = plt.plot(model_dra_list,model_ddec_list,color='lime',label='Solution')
                try:
                    plt.legend((type1,type2,type3,type4,type5[0],type6[0]),('FAST Data Fitted','FAST Data Rejected','NDAC Data Fitted','NDAC Data Rejected','Barycentric motion','Solution'))
                except NameError:
                    try:
                        plt.legend((type1,type3,type4,type5[0],type6[0]),('FAST Data Fitted','NDAC Data Fitted','NDAC Data Rejected','Barycentric motion','Solution'))
                    except NameError:
                        plt.legend((type1,type3,type5[0],type6[0]),('FAST Data Fitted','NDAC Data Fitted','Barycentric motion','Solution'))
                plt.show()
        elif consortium == 'None' or consortium == 'none':
            plt.figure(figsize=(10,10))
            plt.title(f'HIP {HIP}')
            plt.ylabel(r'$\Delta\delta$'+' [mas]')
            plt.xlabel(r'$\Delta\alpha$'+r'$ cos(\delta)$'+' [mas]')
            type5 = plt.plot(barycentric_dra_list,barycentric_ddec_list,color='violet',label='Barycentric motion')
            type6 = plt.plot(model_dra_list,model_ddec_list,color='lime',label='Solution')
            plt.legend((type5[0],type6[0]),('Barycentric motion','Solution'))
            plt.show()
        else:
            print('No such consortium exits. Please check your input consortium.')
    except IndexError:
        print(f'The data of HIP {HIP} cannot be found.')