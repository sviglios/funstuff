# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 14:25:16 2019

@author: konka
"""

import os
import sys
import requests
import time
import argparse
from selenium import webdriver
from bs4 import BeautifulSoup
from Bio import ExPASy
from Bio import SwissProt
import pandas as pd

'''
TO do:
    1. Tests in other PCs
    
human_id = 479
mouse_id = 446
'''


class ExtractDBs:
    """Usage: ExtractDBs(accs, out_folder), where accs is the list of Uniprot accession numbers.
        Extracts peptides from PeptideAtlas Disctinct Observed Peptides table,
        and from ProteomicsDB reference peptide table, for the given accessions.
        Because ProteomicsDB uses client based javascript rendering, the ProteomicsDB
        extraction takes place through a dummy browser, downloading
        the table and processing it. This might not work (altough by iterating 10 times over the page,
        results should be obtained), but Peptide Atlas table should
        be extracted at all times. If the output shows exception for a given accession for 
        the ProteomicsDB database, try to run the script again. Alternatively, the table could be 
        downloaded manually and renamed to
        {acc}refpep.csv, and added to the folder with the rest of the files for further 
        processing."""
        
    def __init__(self, accs, out_folder):
        
        print('\nInitializing')
        self.accs = accs
        self.out_folder = out_folder
        

    def MakeFolder(self):
        
        if self.out_folder not in os.listdir():
            os.mkdir(self.out_folder)
            
        
    def GetPeptidesHTML(self):
        
        options = webdriver.ChromeOptions() 
        options.add_argument("download.default_directory=C:/Downloads")
        driver = webdriver.Chrome(options=options)
        dlpath = r"C:\Users\konka\Downloads/"                  #how do you generalize this without knowing the user?
        destpath = r'C:\Users\konka\Documents\Python Scripts\DBsearch/' + self.out_folder + '/'  #os.getcwd() could work here
        ppflag = 1
        
        for acc in self.accs:
            
            if acc + 'dopAtlas.txt' in os.listdir(self.out_folder):
                print(f'\nFile for {acc} from Peptide Atlas already exists, moving on.')
            else:
                print('Requesting pages for',acc)
                
                url = requests.get("https://db.systemsbiology.net/sbeams/cgi/PeptideAtlas/GetProtein?atlas_build_id=479&protein_name="+acc+"&action=QUERY")
                htmltext = url.text
                
                if not url.ok:
                    print('error', url.status_code)
                
                
                soup = BeautifulSoup(htmltext, 'html.parser')
                
                doptb = soup.find_all(id='getprotein_observedlist_div')
                
                print('Parsing',acc)
                
                if len(doptb) == 1:
                    doptb = doptb[0]
                    
                    doptb_data = doptb.find_all('table')
                    
                    headers = []
                    
                    table1 = doptb_data[0]
                    rt1 = table1.tr
                    
                    ch = 0
                    while ch != 'end':
                        
                        headers.append(rt1.td.text)
                        rt1 = rt1.next_sibling
                        rt1 = rt1.next_sibling
                        if rt1 is None:
                            ch = 'end'
                    
                    table2 = doptb_data[1]
                    data = []
                    
                    rt2 = table2.tr
                    rt2 = rt2.next_sibling
                    rt2 = rt2.next_sibling
                    
                    ch = 0
                    while ch != 'end':
                        
                        rowdata = []
                        celldata = rt2.td
                        
                        ch2 = 0
                        while ch2 != 'end':
                            
                            rowdata.append(celldata.text)
                            
                            celldata = celldata.next_sibling
                            
                            if celldata == '\n':
                                celldata = celldata.next_sibling                
                            if celldata is None:
                                ch2 = 'end'
                                
                        rt2 = rt2.next_sibling
                        
                        if rt2 == '\n':
                            rt2 = rt2.next_sibling            
                        if rt2 is None:
                            ch = 'end'
                        
                        data.append(rowdata)
                            
                else:
                    print("Don't know what to do, more than one tables for", acc)
                    continue
                
                #make peptide accession numbers index, drop their column
                df_dop = pd.DataFrame(data,columns=headers)
                df_dop.index = df_dop.Accession
                df_dop = df_dop.drop('Accession', axis = 1)
                
                #if you want human readable format
                #fh_out = open('Files/'+acc+'dopAtlas.txt','w')
                #fh_out.write(df_dop.to_string(justify='left',index=True))
                #fh_out.close()
                
                df_dop.to_csv(self.out_folder+'/'+acc+'dopAtlas.txt')            
                    
                print(acc, len(data))
                print('Finished PeptideAtlas',acc)
            
            if acc + 'refpep.csv' in os.listdir(self.out_folder):
                print(f'\nFile for {acc} from ProteomicsDB already exists, moving on.')
            else:
                try:
                    handle = ExPASy.get_sprot_raw(acc)
                    record = SwissProt.read(handle)
                except:
                    print('Not found')
                    
                for ref in record.cross_references:
                    
                    ch = 0
                    if ref[0] == 'ProteomicsDB':
                        print('\nProteomicsDB', acc, ref[1])
                        protid = ref[1]
                        ch = 1
                    
                    if ch == 1:
                        
                        end = 0
                        for i in range(10):
                            driver.get('https://www.proteomicsdb.org/proteomicsdb/#human/proteinDetails/' + protid + '/referencePeptides')
                            button = driver.find_element_by_xpath('/html/body/div[1]/div/div[3]/article/div/table/tbody/tr[4]/td/div/div/div/div[2]/div/table/tbody/tr[2]/td/div/div[1]/div/div/div/button/span[2]')
                            print('Iteration', i+1)
                            try:
                                button.click()
                                time.sleep(3)
                                end = 1
                                break
                            except:
                                time.sleep(1)
                                continue
                        
                        if end == 0:
                            print('EXCEPTION FOR', acc, protid, '\n\nRun the script again, so the extraction can be attempted again. No post-processing will occur. If you want to activate it anyway, pass arg for postprocessing\n\n')
                            ppflag = 0
                            continue
                        
                        print('Downloaded peptides for', acc)
                        
                        if acc + 'refpep.csv' not in os.listdir(self.out_folder):
                            os.rename(dlpath+'peptides.csv', destpath + acc + 'refpep.csv')
                        else:
                            print(f'File for {acc} already exists, moving on.')
                            os.remove(r"C:\Users\konka\Downloads/" + 'peptides.csv')
                        
                        print('Finished ProteomicsDB', acc, '\n\n')
                        break
            
        time.sleep(2)
        driver.quit()
        
        return ppflag
    
    
    def UniprotFind(acc):
        
        try:
            handle = ExPASy.get_sprot_raw(acc)
            record = SwissProt.read(handle)
            print('\nUniprot data extracted for', acc)
            #sequence = record.sequence    
        except:
            return None, None, None
        
        flag = 0
        for feat in record.features:
            if feat.type == 'SIGNAL':
                end = list(feat.location)[-1]
                flag = 1
            if feat.type == 'PROPEP':
                end = list(feat.location)[-1]
                flag = 1
                break
        
        if flag == 0:
            return record.entry_name, 0, record.sequence_length
        
        return record.entry_name, end, record.sequence_length
    
    
    def PeptideFilter(pep):
        
        pos = len(pep) - 1
        fl = 0
        
        for a in range(len(pep)):
            if len(pep) < 8 or len(pep) > 21:
                fl = 1
            if pep[a] == 'R' or pep[a] == 'K':
                if pos != a:
                    fl = 2
            if pep[a] == 'M':
                fl = 3
            if a == pos:
                if pep[a] != 'K' and pep[a] != 'R':
                    fl = 4
        
        return fl
        
    
    def ProcessTables(self):
        
        dfs = []
        
        for acc in self.accs:
            
            name, pro_end, length = ExtractDBs.UniprotFind(acc)
            if name == None:
                print('\n', acc, 'NOT FOUND in Uniprot, skiping!\n')
                continue
            else:
                print(name, pro_end, length)
            
            df_atlas = pd.read_csv(self.out_folder + '/' + acc + 'dopAtlas.txt')
            for i in range(len(df_atlas)):
                df_atlas.loc[df_atlas.index[i],'ESS'] = df_atlas.loc[df_atlas.index[i],'ESS'].split(' ')[0]    
            df_atlas = df_atlas.sort_values('ESS', ascending = False)
            
            try:
                df_pdb = pd.read_csv(self.out_folder + '/' + acc + 'refpep.csv', sep=';')
                df_pdb = df_pdb.sort_values('ANDROMEDA_SCORE', ascending = False)
                
                df_com = pd.merge(df_atlas, df_pdb, how='inner', left_on='Sequence', right_on ='PLAIN_SEQUENCE')
                df_com = df_com.sort_values(['MISSED_CLEAVAGES','ANDROMEDA_SCORE','ESS'], ascending = [True,False,False])
                
            except:
                print('\nProteomicsDB file not available for',acc)
                df_com = df_atlas
            
            df_com = df_com.loc[~df_com.Accession.duplicated(keep='first')]
            
            dec = 0
            for i in range(len(df_com)):
                i  -= dec
                ind = df_com.index[i]
                pep = df_com.loc[ind,'Sequence']
                
                flag = ExtractDBs.PeptideFilter(pep)
                if flag in [1,2,4]:
                    df_com = df_com.drop(ind)
                    print('Dropped',pep,'flag =', flag)
                    dec += 1
                elif flag == 3:
                    if len(df_com) > 5:
                        df_com = df_com.drop(ind)
                        print('Dropped',pep,'flag = ', 3)
                        dec += 1
            
            print(f'\nFinal number of peptides for {acc} is {len(df_com)}, picking 15 or full dataframe if less')            
            df_sel = df_com.iloc[:15,:]
            
            ind = [acc+'-'+name for i in range(len(df_sel))]
            df_sel = df_sel.set_index([ind,'Accession'])    
            dfs.append(df_sel)
                
            print('\nFinished', acc)
            
        df_final = pd.concat(dfs, sort = True)
        df_final.to_excel(self.out_folder + '/' + 'FinalPeps.xlsx',index=True)
        
        fh_out = open(self.out_folder + '/' + 'FinalPeps.txt','w')
        fh_out.write(df_final.to_string(justify='left',index=True))
        fh_out.close()
        
        return dfs
    
    
    def MouseHomolog(self, dfs):
        
        print('\nFinding mouse homologs')
        ind = 0
        new_dfs = []
        
        for acc in self.accs:
            
            try:
                handle = ExPASy.get_sprot_raw(acc)
                record = SwissProt.read(handle)
                name = record.entry_name    
            except:
                print('\nNo entry for', acc, ',continuing')
                ind += 1
                continue
            
            try:
                mname = name.split('_')[0] + '_MOUSE'
                mhandle = ExPASy.get_sprot_raw(mname)
                mrecord = SwissProt.read(mhandle)
                mseq = mrecord.sequence
                print(f'\nFound mouse homolog for {name}: {mname}')
            except:
                print(f'\nNo mouse gene entry for {acc}-{name}, continuing')
                ind += 1
                continue
            
            df = dfs[ind]
            mcol = []
            
            for row in range(len(df)):
                pepseq = df.Sequence[df.index[row]]
                print(pepseq)
                if str(pepseq) in mseq:
                    mcol.append('True')
                else:
                    mcol.append('False')
            
            df['Mouse'] = mcol
            new_dfs.append(df)
            ind += 1
            
        df_final = pd.concat(new_dfs, sort = True)
        df_final.to_excel(self.out_folder + '/' + 'MouseHomologPeptides.xlsx',index=True)


    def search_EBI(self):
        
        if 'EBI' not in os.listdir(self.out_folder):
            os.mkdir(self.out_folder + '/EBI')
            
        for acc in self.accs:
            #requestURL = "https://www.ebi.ac.uk/proteins/api/proteomics?offset=0&size=100&accesion=P14780&taxid=9606&datasource=PeptideAtlas&unique=True"    
            requestURL = "https://www.ebi.ac.uk/proteins/api/proteomics/" +acc
            
            try:
                r = requests.get(requestURL, headers={ "Accept" : "application/json"})
            except:
                print(r.status_code)
                print(acc)
                continue
            
            if not r.ok:
                #r.raise_for_status()
                print(r.status_code)
                print(acc)
                continue
            
            res = r.json()
            r.close()
            
            print('\nEBI request ok, parsing\n')
            
            data = res['features']
            
            df_eb = pd.DataFrame.from_dict(data)
            dbcol = []

            for i in range(len(df_eb)):
                dbna = []
                for d in range(len(df_eb.loc[i]['evidences'])):
                    dbna.append(df_eb.loc[i]['evidences'][d]['source']['name'])
                dbcol.append(','.join(dbna))
      
            df_out3 = df_eb.loc[:, ['begin', 'end', 'peptide', 'unique']]
            df_out3.loc[:, 'database'] = dbcol
            
            print(f"\nWriting file for {res['entryName']}\n")
            
            filename = self.out_folder + '/EBI/' + acc + 'ebi.xlsx'
            df_out3.to_excel(filename,index=False)
    
    
    def search_protDB(self):
        
        if 'ProtDB' not in os.listdir(self.out_folder):
            os.mkdir(self.out_folder + '/ProtDB')
        
        for acc in self.accs:
            resp = requests.get("https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinpeptideresult.xsodata/InputParams(PROTEINFILTER='"+acc+"')/Results?$select=ENTRY_NAME,PROTEIN_NAME,UNIQUE_IDENTIFIER,TAXCODE,CHROMOSOME_NAME,GENE_NAME,STRAND,PEPTIDE_SEQUENCE,PEPTIDE_MASS,START_POSITION,END_POSITION,SCORE,RANK,Q_VALUE,PEP,SEARCH_ENGINE,ISUNIQUE,ISUNIQUE_PROTEIN,PROJECT_NAME,PROJECT_DESCRIPTION,EXPERIMENT_NAME,EXPERIMENT_ID,EXPERIMENT_DESCRIPTION,PUBMEDID,SCOPE&$filter=PEPTIDE_MASS%20gt%201000%20&$format=json")
            
            if resp.status_code == 200:
                print('\nRecieved response, ProtDB request successful for', acc)
            else:
                print('\nSomething went wrong, try again')
                sys.exit()
            
            print('\nParsing')
            
            result = resp.json()
            resp.close()
            
            res_dic = result['d']['results']
            
            df = pd.DataFrame.from_dict(res_dic)
            
            try:
                df.loc[:, "SCORE"] = df.loc[:, "SCORE"].astype(float)
            except KeyError:
                print(f'No proteomics file for {acc}')
                print(df)
                continue
            
            df_un = df[df['ISUNIQUE_PROTEIN'] == 1 ]
            df_un = df_un.sort_values("SCORE", ascending=False)
            
            df_un.loc[df_un['SEARCH_ENGINE'] == 2, 'SEARCH_ENGINE'] = 'MASCOT'
            df_un.loc[df_un['SEARCH_ENGINE'] == 1, 'SEARCH_ENGINE'] = 'ANDROMEDA'
            
            df_out = df_un[['EXPERIMENT_ID','PEPTIDE_SEQUENCE','PEP','Q_VALUE','RANK','SCORE','SEARCH_ENGINE']]
            
            filename = self.out_folder + '/ProtDB/' + acc +'protDB.xlsx'
            df_out.to_excel(filename,index=False)
            print('\nWrote file:',filename)
            #time.sleep(3)
            
            resp = requests.get("https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinproteotypicity.xsodata/InputParams(PROTEINFILTER='"+acc+"',LABEL_TYPES='SILAC')/Results?$select=UNIQUE_IDENTIFIER,RANK_ORDER,PEPTIDE_ID,PSMS,OCCURRENCE,PROTEOTYPICITY,CUM_PROTEOTYPICITY,GENE_COUNT,IDENTIFIER_COUNT,SEQUENCE,UNIQUENESS&$format=json")    
            if resp.status_code == 200: #a better way could also be {if resp.ok:}
                print('\nRecieved response, proteotypicity request successful for', acc)
            else:
                print('\nSomething went wrong, try again')
                print(resp.status_code)
                sys.exit()
            
            print('\nParsing')
            
            result = resp.json()
            resp.close()
            
            res_dic = result['d']['results']
            
            df = pd.DataFrame.from_dict(res_dic)
            
            try:
                df["CUM_PROTEOTYPICITY"] = df["CUM_PROTEOTYPICITY"].astype(float)    
            except KeyError:
                print(f'No proteotypicity file for {acc}')
                print(df)
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
                
                df_pepsel = df_pep.loc[:, ['SCORE','INTENSITY','Q_VALUE']]

                cols = []
                dat_l = []
                for i in df_pepsel.columns:
                    cols.append(i)
                    df_pepsel.loc[:, i] = df_pepsel.loc[:, i].astype(float)
                    dat_l.append(df_pepsel[i].mean())
   
                ser_l.append(dat_l)
            
            df_out3 = pd.DataFrame(ser_l,columns = cols)

            for i in df_out3.columns:
                df_out2.loc[:, i] = df_out3.loc[:,i]
            
            filename =  self.out_folder + '/ProtDB/' + acc + 'refpepprotDB.xlsx'
            df_out2.to_excel(filename,index=False)
            print('\nWrote file:',filename)


def ReadFile(filename):
    
    accs = []
    fh = open(filename,'r')
    
    for line in fh:
        
        line = line.strip()
        accs.append(line)
    
    fh.close()
    
    return accs


########################## TEST ##########################

# =============================================================================
# obex = ExtractDBs(['P01033','P16035','P35625','Q99727'], 'PeptidesClass')
# 
# obex.MakeFolder()
# 
# ppflag = obex.GetPeptidesHTML()
# 
# if ppflag == 1:
#     dfs = obex.ProcessTables()
#     obex.MouseHomolog(dfs)
#     print('\n\nFinished, closing down.\n')
# else:
#     print('\nPost processing omitted, shutting down.\n')
#     
# obex.search_EBI()
# obex.search_protDB()
# =============================================================================

        
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description = 'This class searches for \
                                     proteotypic for the proteins provided (either\
                                     through a file handle or through the command. Two\
                                     different ways of extracting the peptides\
                                     are employed, web scraping and api calls.\
                                     Web scraping is the default option, add -s\
                                     to also query the APIs of EBI and ProteomicsDB.\
                                     Call the script with -m to also search for\
                                     peptides that are homolog in mouse.', 
                                     epilog = 'Flags -s and -m do not need arguments,\
                                     but the value True or False can be supplied, although \
                                     it is not necessary.', 
                                     argument_default = None,
                                     add_help = True, 
                                     formatter_class=argparse.HelpFormatter, 
                                     prefix_chars = '-') #optional args prefix
    #None is already the default value of the argument_default argument HA-HA
    
    parser.add_argument('-f', '--filename', default = None, help = 'A file \
                        input with each accession number in different line')
    
    parser.add_argument('-o', '--outfolder', default = 'PeptidesDBs',
                        help = 'Destination folder for files', nargs = '?')
    
    parser.add_argument('-m', '--mouse', default = False, help = 'Check \
                        if they peptides are present in mouse homolog proteins.\
                        Default is False', type = bool, const = True, nargs = '?')
    
    parser.add_argument('-s', '--search', default = False, help = 'Search in \
                        EBI and ProteomicsDB databases through the API. \
                        Default is False', type = bool, const = True, nargs = '?')
    
    parser.add_argument('-a', '--accs', default = None, help = 'Add accession\
                        number through the command line instead of a file.', 
                        type = str, nargs='+')
    
    args = parser.parse_args()
    
    if args.filename and args.accs:
        print('Both a filename and accessions from the command line were supplied. Reading from file.')        
        accs = ReadFile(args.filename)
        
    elif args.filename:
        accs = ReadFile(args.filename)
    
    elif args.accs:
        accs = args.accs
    
    else:
        print('No accessions supplied. Use -f to provide a file or -a to provide raw accessions from cmd. Exiting..')
        time.sleep(4)
        sys.exit(4)
    
    obex = ExtractDBs(accs, args.outfolder)
    obex.MakeFolder()
    ppflag = obex.GetPeptidesHTML()
    
    if ppflag == 1:
        dfs = obex.ProcessTables()
        
        if args.mouse:
            obex.MouseHomolog(dfs)
        
        if args.search:
            obex.search_EBI()
            obex.search_protDB()
        
        print('\n\nFinished, closing down.\n')
    
    else:
        print('\nPost processing omitted, shutting down.\n')
    
    time.sleep(2)


########################## END TEST ##########################