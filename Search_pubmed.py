# -*- coding: utf-8 -*-
"""
Created on Mon Jul  1 12:38:14 2019
@author: Kostas
"""

#Biopython can be installed with:
'''pip install biopython'''

import sys
import time
import pandas as pd
import datetime
from Bio import Entrez
from Bio import Medline

print('This script retrieves pubmed IDs from the term of interest\n')

#input vars
Entrez.email = input('Email (NCBI requirements): ') # Always tell NCBI who you are
database = input('NCBI Database to search: ')
term = input('Term to search: ')
maxres = input('Max number of returned results: ')
mindate = input('Date lower limit (press enter for default): ')
maxdate = input('Date upper limit (press enter for default): ')

# =============================================================================
# #testing
# Entrez.email = 'example@email.com'
# database = 'pubmed'
# term = 'heart disease'
# maxres = 25
# mindate = 2001
# maxdate = 2005
# =============================================================================

if mindate == '' or mindate == 0:
    mindate = None
if maxdate == '' or maxdate == 0:
    maxdate = None
    
#search terms for entries adn get their ids, raise error if database input not correct
try:
    handle = Entrez.esearch(db=database, term=term, retmax = maxres, mindate=mindate, maxdate=maxdate, idtype='acc')
    result = Entrez.read(handle)
    handle.close()
except RuntimeError:
    print('\nDatabase name not valid, please try again\n')
    print('Exiting in 5 seconds')
    time.sleep(5)
    sys.exit(4)
    
idlist = result["IdList"]

filenametxt = '_'.join(term.split()) + '.txt'

if database != 'pubmed':
    fh = open(filenametxt,'w')
    for i in idlist:
        fh.write(str(i)+'\n')
    fh.close()

else:    
    print(f'\nExtracted {len(idlist)} ids. Database is pubmed, requesting records.')
    
    print('Getting entries records..')
    
    handle = Entrez.efetch(db='pubmed', id=idlist, rettype='medline', retmode='text')
    records = Medline.parse(handle)
    records = list(records)
     
    filenamexl = '_'.join(term.split()) + '.xls'
    fh = open(filenametxt,'w')
    
    cols = ['Pubmed ID','Date','Type','Journal','Citations','TITLE','Authors'] 
    data = []
    
    c = 0
    v = 10
    #write data
    for rec in records:
        
        try:
            date = datetime.datetime.strptime(rec['PHST'][-1].split()[0], '%Y/%m/%d').strftime('%d/%m/%y')
        except KeyError:
            print('Search yielded no results. Please change your query and try again')
            fh.close()
            print('Exiting in 5 seconds')
            time.sleep(5)
            sys.exit(3)
        
        citations = Entrez.read(Entrez.elink(dbfrom="pubmed", db="pmc", LinkName="pubmed_pmc_refs", id = rec['PMID']))
        
        #catch error
        if len(citations[0]["LinkSetDb"]) == 0:
            pmc_ids = []
        else:
            pmc_ids = [link["Id"] for link in citations[0]["LinkSetDb"][0]["Link"]]
        
        if 'AU' not in rec:
            rec['AU'] = 'No authors'
            
        data.append([rec['PMID'], 
                     date,
                     '/'.join(rec['PT']),
                     rec['JT'],
                     len(pmc_ids),
                     rec['TI'], 
                     ' ,'.join(rec['AU'])])
        
        perc = 100*(c/len(records))
        if perc > v:
            print(f'{perc}% done...')
            v += 10
        c += 1
    
    print('Writing file...')
    
    #make df to write
    pd.set_option("display.max_colwidth", 150)
    df = pd.DataFrame(data,columns=cols)
    df.index = [i+1 for i in list(df.index)]
    
    #sort by citations
    df = df.sort_values("Citations", ascending=False)
    
    fh.write(df.to_string(justify='justify-all'))
    
    fh.close()
    
    df.to_excel(filenamexl,index=False)

#finish statements
print("\nDone\n")
print("Text file is named: " + filenametxt)

if database == 'pubmed':
    print("Excel file is named: " + filenamexl)

print(f"\nClosing in 10 seconds")
time.sleep(10)

#finish with countdown
# =============================================================================
# t = 10
# for i in range(10):
#     print(f"\nClosing in {t} seconds")
#     t -= 1
#     time.sleep(1)
# =============================================================================