# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 09:58:56 2019

@author: Kostas
"""

import requests
import time
import pandas as pd
import sys

accs = ['P14780','P08253']

accs = []
fh = open('prot.txt','r')
for line in fh:
    line = line.strip()
    accs.append(line)

for acc in accs:
    resp = requests.get("https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinpeptideresult.xsodata/InputParams(PROTEINFILTER='"+acc+"')/Results?$select=ENTRY_NAME,PROTEIN_NAME,UNIQUE_IDENTIFIER,TAXCODE,CHROMOSOME_NAME,GENE_NAME,STRAND,PEPTIDE_SEQUENCE,PEPTIDE_MASS,START_POSITION,END_POSITION,SCORE,RANK,Q_VALUE,PEP,SEARCH_ENGINE,ISUNIQUE,ISUNIQUE_PROTEIN,PROJECT_NAME,PROJECT_DESCRIPTION,EXPERIMENT_NAME,EXPERIMENT_ID,EXPERIMENT_DESCRIPTION,PUBMEDID,SCOPE&$filter=PEPTIDE_MASS%20gt%201000%20&$format=json")
    
    if resp.status_code == 200:
        print('\nRecieved response, request successful')
    else:
        print('\nSomething went wrong, try again')
        time.sleep(3)
        sys.exit()
    
    print('\nParsing')
    
    result = resp.json()
    resp.close()
    
    res_dic = result['d']['results']
    
    df = pd.DataFrame.from_dict(res_dic)
    df["SCORE"] = df["SCORE"].astype(float)
    
    df_un = df[df['ISUNIQUE_PROTEIN'] == 1 ]
    df_un = df_un.sort_values("SCORE", ascending=False)
    
    df_un.loc[df_un['SEARCH_ENGINE'] == 2, 'SEARCH_ENGINE'] = 'MASCOT'
    df_un.loc[df_un['SEARCH_ENGINE'] == 1, 'SEARCH_ENGINE'] = 'ANDROMEDA'
    
    df_out = df_un[['EXPERIMENT_ID','PEPTIDE_SEQUENCE','PEP','Q_VALUE','RANK','SCORE','SEARCH_ENGINE']]
    
    filename = acc +'protDB.xlsx'
    df_out.to_excel(filename,index=False)
    print('\nWrote file:',filename)
    #time.sleep(3)
    
    resp = requests.get("https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinproteotypicity.xsodata/InputParams(PROTEINFILTER='"+acc+"',LABEL_TYPES='SILAC')/Results?$select=UNIQUE_IDENTIFIER,RANK_ORDER,PEPTIDE_ID,PSMS,OCCURRENCE,PROTEOTYPICITY,CUM_PROTEOTYPICITY,GENE_COUNT,IDENTIFIER_COUNT,SEQUENCE,UNIQUENESS&$format=json")    
    if resp.status_code == 200: #a better way could also be {if resp.ok:}
        print('\nRecieved response, request successful')
    else:
        print('\nSomething went wrong, try again')
        print(resp.status_code)
        #time.sleep(3)
        sys.exit()
    
    print('\nParsing')
    
    result = resp.json()
    resp.close()
    
    res_dic = result['d']['results']
    
    df = pd.DataFrame.from_dict(res_dic)
    
    try:
        df["CUM_PROTEOTYPICITY"] = df["CUM_PROTEOTYPICITY"].astype(float)    
    except KeyError:
        print(f'No file for {acc}')
        continue

    df = df.sort_values("CUM_PROTEOTYPICITY", ascending=False)
    
    df_out2 = df[['CUM_PROTEOTYPICITY', 'OCCURRENCE', 'PEPTIDE_ID', 'PROTEOTYPICITY', 'PSMS', 'RANK_ORDER', 'SEQUENCE', 'UNIQUENESS']]
    
    df_out2.index = range(len(df_out2.index))
    
    ser_l = []
    for seq in df_out2['SEQUENCE']:
        resp = requests.get("https://www.proteomicsdb.org/proteomicsdb/logic/api/peptidesearch.xsodata/InputParams(PEPTIDE_SEQUENCE='" + seq +"',Q_VALUE_CUTOFF=0.01)/Results?$select=IDENTIFICATION_ID,SEQUENCE,SCORE,PEPTIDE_ID,Q_VALUE,PROTEASE_NAME,SEARCH_ENGINE_NAME,QUANTIFICATION_METHOD_ID,INTENSITY,QUANTIFICATION_TYPE&$format=json")
        if resp.ok:
            print(f'Peptide {seq} successful')
        else:
            print(f'Peptide search {seq} failed')
        
        res = resp.json()
        resp.close()
        
        res_dic = res['d']['results']
        
        df_pep = pd.DataFrame.from_dict(res_dic)
        
        df_pepsel = df_pep[['SCORE','INTENSITY','Q_VALUE']]
        
        cols = []
        dat_l = []
        for i in df_pepsel.columns:
            cols.append(i)
            df_pepsel[i] = df_pepsel[i].astype(float)
            dat_l.append(df_pepsel[i].mean())
            
        ser_l.append(dat_l)
    
    df_out3 = pd.DataFrame(ser_l,columns = cols)
    
    for i in df_out3.columns:
        df_out2[i] = df_out3[i]
    
    filename = acc + 'refpepprotDB.xlsx'
    df_out2.to_excel(filename,index=False)
    print('\nWrote file:',filename)
    #time.sleep(3)
