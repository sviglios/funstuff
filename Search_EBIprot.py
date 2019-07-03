# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 15:28:04 2019

@author: Kostas
"""

import requests, sys, pandas as pd

accs = ['P14780','P08253']

accs = []
fh = open('prot.txt','r')
for line in fh:
    line = line.strip()
    accs.append(line)

for acc in accs:
    #requestURL = "https://www.ebi.ac.uk/proteins/api/proteomics?offset=0&size=100&accesion=P14780&taxid=9606&datasource=PeptideAtlas&unique=True"    
    requestURL = "https://www.ebi.ac.uk/proteins/api/proteomics/" +acc
    r = requests.get(requestURL, headers={ "Accept" : "application/json"})
    
    if not r.ok:
        r.raise_for_status()
        sys.exit()
    
    res = r.json()
    r.close()
    
    print('\nRequest ok, parsing\n')
    
    data = res['features']
    datad = {}
    
    df_eb = pd.DataFrame.from_dict(data)
    dbcol = []
    
    for i in range(len(df_eb)):
        dbna = []
        for d in range(len(df_eb.loc[i]['evidences'])):
            dbna.append(df_eb.loc[i]['evidences'][d]['source']['name'])
        dbcol.append(','.join(dbna))
            
    df_out3 = df_eb[['begin', 'end', 'peptide', 'unique']]
    df_out3['database'] = dbcol
    
    print(f"\nWriting file for {res['entryName']}\n")
    
    filename = acc +'ebi.xlsx'
    df_out3.to_excel(filename,index=False)