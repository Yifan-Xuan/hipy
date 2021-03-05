import requests

def get_data(id,flag):
  if flag == 1:
    url="https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?hipId="+str(id)
    r = requests.get(url)
    with open('/home/yifan_xuan/catalogue/'+str(id)+'_catalog.txt','wb') as f:
        f.write(r.content)
  elif flag == 2:
    url="https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?noLinks=1&tabular=1&hipiId="+str(id)
    r = requests.get(url)
    with open('/home/yifan_xuan/intermediate_data/'+str(id)+'_Id.txt','wb') as f:
        f.write(r.content)
  elif flag == 3:
    url="https://hipparcos-tools.cosmos.esa.int/cgi-bin/HIPcatalogueSearch.pl?hipepId="+str(id)
    r = requests.get(url)
    with open('/home/yifan_xuan/epoch_photometry_data/'+str(id)+'_epd.txt','wb') as f:
        f.write(r.content)
  return

