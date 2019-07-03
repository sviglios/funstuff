# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 10:06:03 2019

@author: Kostas
"""

import pandas as pd

accs = []
fh = open('prot.txt','r')
for line in fh:
    line = line.strip()
    accs.append(line)

fh_out = open('Final_peptides.txt','w')

count = 1    
for acc in accs:
    
    fh_out.write(f'{count}. Peptides for {acc}')
    eb = 1
    pdb = 1
    
    try:
        df_ebi = pd.read_excel('Files/Sim/EBI_DB/'+acc+'ebi.xlsx')
        df_ebi = df_ebi.loc[~df_ebi.peptide.duplicated(keep='first')]
    except:
        print(f'No file for {acc} in EBI')
        eb = 0
    
    try:
        df_pdb = pd.read_excel('Files/Sim/RefProtDB/'+acc+'refpepprotDB.xlsx')
        df_pdb = df_pdb.loc[~df_pdb.SEQUENCE.duplicated(keep='first')]
    except:
        print(f'No file for {acc} in protDB')
        pdb = 0
    
    if eb == 0 and pdb == 0:
        
        print(f'\nCould not find peptides for the accession {acc}')
        fh_out.write(f'Could not find peptides for the accession {acc}')
        
    elif eb == 0 and pdb == 1:
        
        print(f'That is weird, {acc} is present in proteomics DB but not in EBI. Skipping this protein.')
        fh_out.write(f'That is weird, {acc} is present in proteomics DB but not in EBI. Skipping this protein.')
    
    elif eb == 1 and pdb == 0:
        
        print(f'Cound not find peptides in proteomics DB for {acc}. Picking from EBI only')
        fh_out.write(f'\nCound not find peptides in proteomics DB for {acc}. Picking from EBI only')
        
        df_ebisel = df_ebi[df_ebi['unique'] == True]
        df_ebisel['db number'] = [len(sp.split(',')) for sp in df_ebisel['database']]
        df_ebisel2 = df_ebisel[df_ebisel['db number'] == 4]
        cutoff = 3
        
        while len(df_ebisel2) < 3:
            print(f'Less than 3 peptides in EBI, accession {acc}')
            df_ebisel2 = df_ebisel[df_ebisel['db number'] >= cutoff]
            print(f'{len(df_ebisel2)} for {acc} including peptides with {cutoff} databases')
            cutoff -= 1
        else:
            print(f'{acc}\t{len(df_ebisel2)}')
        
        maxdf = max(df_ebisel2['end'])
        mindf = min(df_ebisel2['begin'])
        
        ch = 0
        c = 0         
        while len(df_ebisel2) > 3:            
            if ch == 0:
                df_ebisel2 = df_ebisel2.drop(df_ebisel2[df_ebisel2['end'] == maxdf].index)
                print('Removed end')
                ch = 1
            elif ch == 1:
                df_ebisel2 = df_ebisel2.drop(df_ebisel2[df_ebisel2['begin'] < 20].index)
                print('Removed start')
                ch = 2
            elif ch == 2:
                if c < 10:
                    df_ebisel2 = df_ebisel2.sort_values(['begin'])
                    print(min(df_ebisel2['begin']))
                    df_ebisel2 = df_ebisel2.drop(df_ebisel2.index[0])
                    print('Removed instance, new begin is',min(df_ebisel2['begin']))
                    c += 1
                else:
                    df_ebisel2 = df_ebisel2.sample(3)
        
        print(acc,' ',len(df_ebisel2))
        fh_out.write('\n\n')
        fh_out.write(df_ebisel2.to_string(justify='left',index=False))
        fh_out.write('\n===================================================\n\n')
        
    elif eb == 1 and pdb == 1:

        df_ebisel = df_ebi[df_ebi['unique'] == True]
        df_ebisel['db number'] = [len(sp.split(',')) for sp in df_ebisel['database']]
        df_ebisel2 = df_ebisel[df_ebisel['db number'] == 4]
        cutoff = 3        
        
        while len(df_ebisel2) < 3:
            print(f'Less than 3 peptides in EBI, accession {acc}')
            df_ebisel2 = df_ebisel[df_ebisel['db number'] >= cutoff]
            print(f'{len(df_ebisel2)} for {acc} including peptides with {cutoff} databases')
            cutoff -= 1
        print(f'{acc}\t{len(df_ebisel2)}')
        
        seq1 = set(df_ebisel2['peptide'])
        seq2 = set(df_pdb['SEQUENCE'])
        
        new_seq = seq1.intersection(seq2)
            
        if len(new_seq) < 3:            
            
            maxdf = max(df_ebisel2['end'])
            mindf = min(df_ebisel2['begin'])
            
            ch = 0
            c = 0         
            while len(df_ebisel2) > 3:            
                if ch == 0:
                    df_ebisel2 = df_ebisel2.drop(df_ebisel2[df_ebisel2['end'] == maxdf].index)
                    print('Removed end')
                    ch = 1
                elif ch == 1:
                    df_ebisel2 = df_ebisel2.drop(df_ebisel2[df_ebisel2['begin'] == mindf].index)
                    print('Removed start')
                    ch = 2
                elif ch == 2:
                    if c < 10:
                        df_ebisel2 = df_ebisel2.sort_values(['begin'])
                        print(min(df_ebisel2['begin']))
                        df_ebisel2 = df_ebisel2.drop(df_ebisel2.index[0])
                        print('Removed instance, new begin is',min(df_ebisel2['begin']))
                        c += 1
                    else:
                        df_ebisel2 = df_ebisel2.sample(3)
                
                        print(acc,' ',len(df_ebisel2))
            
            fh_out.write('\n\n')
            fh_out.write(df_ebisel2.to_string(justify='left',index=False))
            fh_out.write('\n===================================================\n\n')
        
        else:            
            df_both = df_ebisel2[df_ebisel2.peptide.isin(list(new_seq))]
            df_both['proteotypic'] = list(df_pdb[df_pdb.SEQUENCE.isin(list(new_seq))]['CUM_PROTEOTYPICITY'])
            df_both['uniqueness'] = list(df_pdb[df_pdb.SEQUENCE.isin(list(new_seq))]['UNIQUENESS'])
            maxdf = max(df_both['end'])
            mindf = min(df_both['begin'])
            
            ch = 0
            while len(df_both) > 3:
                if ch == 0:
                    df_both = df_both.drop(df_both[df_both['end'] == maxdf].index)
                    print('Removed end')
                    ch = 1
                elif ch == 1:
                    df_both = df_both.drop(df_both[df_both['begin'] < 20].index)
                    print('Removed start')
                    ch = 2
                elif ch == 2:
                    seqs = list(df_both['peptide'])
                    df_pdb = df_pdb[df_pdb.SEQUENCE.isin(seqs)]
                    df_pdb = df_pdb.sort_values('PROTEOTYPICITY',ascending=False)
                    sel = list(df_pdb['SEQUENCE'][:3])
                    df_both = df_both[df_both.peptide.isin(sel)]
                
            print(acc,' ',len(df_both))
            fh_out.write('\n\n')
            fh_out.write(df_both.to_string(justify='left',index=False))
            fh_out.write('\n===================================================\n\n')
    count += 1
                    