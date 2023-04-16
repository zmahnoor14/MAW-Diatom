#!/usr/bin/env python
# coding: utf-8

# In[276]:


## load libraries


# In[337]:


import glob
import json
import os
import re
import time
import wget
import urllib.parse
import argparse
import plotly.express as px

import numpy as np
import pandas as pd
import pubchempy as pcp


from pybatchclassyfire import *
from pandas import json_normalize
from rdkit import Chem
from rdkit import DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFMCS
from rdkit.Chem import PandasTools
def isNaN(string):
    return string != string


# In[338]:


## Functions


# In[409]:


# naming is any string used to name or put an id to the sunburst as many sunbursts will be created
def sunburst(input_dir, input_csv, naming):
    
    cl = pd.read_csv(input_csv)
    class_data = cl[['superclass', 'class', 'subclass']]
    spclass = list(class_data['superclass']) # all superclasses
    uniq_spclass = list(np.unique(list(class_data['superclass']))) # only unique super classes
    uniq_spc = [s for s in uniq_spclass if 'nan' not in s ] # only unique super classes with no NA values
    print(len(uniq_spclass))
    clss = list(class_data['class'])
    uniq_class = list(np.unique(list(class_data['class'])))
    uniq_c = [s for s in uniq_class if 'nan' not in s ]
    len(uniq_class)
    sbclass = list(class_data['subclass'])
    uniq_sbclass = list(np.unique(list(class_data['subclass'])))
    uniq_sbc = [s for s in uniq_sbclass if 'nan' not in s ]
    len(uniq_sbclass)

    #all characters
    Names = ['Organic Compounds'] + uniq_spclass+uniq_class+uniq_sbclass

    df = pd.DataFrame(Names)
    df['values'] = ''
    df['parents'] = ''

    df = df.rename(columns={0: 'characters'})
    
    if "nan" in np.unique(df["characters"]):
    
        #for i, row in df.iterrows():
        for i, row in df[0:len(df)-2].iterrows():
            if 'Organic Compounds' in df['characters'][i]:
                df.loc[i, 'values'] = 0
                df.loc[i, 'parents'] = ''

            elif df['characters'][i] in uniq_spclass:

                df.loc[i, 'values'] = spclass.count(df['characters'][i])
                df.loc[i, 'parents'] = 'Organic Compounds'

            elif df['characters'][i] in uniq_class:

                df.loc[i, 'values'] = clss.count(df['characters'][i])
                df.loc[i, 'parents'] = 'Organic Compounds'

                df.loc[i, 'values'] = clss.count(df['characters'][i])
                clsp = class_data['superclass'][class_data[class_data['class'] == df['characters'][i]].index.tolist()[0]]
                df.loc[i, 'parents'] = clsp


            elif df['characters'][i] in uniq_sbclass:
                df.loc[i, 'values'] = sbclass.count(df['characters'][i])
                sbclsp = class_data['class'][class_data[class_data['subclass'] == df['characters'][i]].index.tolist()[0]]
                df.loc[i, 'parents'] = sbclsp
    else:
        for i, row in df.iterrows():
            if 'Organic Compounds' in df['characters'][i]:
                df.loc[i, 'values'] = 0
                df.loc[i, 'parents'] = ''

            elif df['characters'][i] in uniq_spclass:

                df.loc[i, 'values'] = spclass.count(df['characters'][i])
                df.loc[i, 'parents'] = 'Organic Compounds'

            elif df['characters'][i] in uniq_class:

                df.loc[i, 'values'] = clss.count(df['characters'][i])
                df.loc[i, 'parents'] = 'Organic Compounds'

                df.loc[i, 'values'] = clss.count(df['characters'][i])
                clsp = class_data['superclass'][class_data[class_data['class'] == df['characters'][i]].index.tolist()[0]]
                df.loc[i, 'parents'] = clsp


            elif df['characters'][i] in uniq_sbclass:
                df.loc[i, 'values'] = sbclass.count(df['characters'][i])
                sbclsp = class_data['class'][class_data[class_data['subclass'] == df['characters'][i]].index.tolist()[0]]
                df.loc[i, 'parents'] = sbclsp
    data = dict(character = df['characters'], parents = df['parents'], values = df['values'])
    fig = px.sunburst(
        data,
        names='character',
        parents='parents',
        values='values',
        
    )
    fig.update_layout(margin = dict(t=0, l=0, r=0, b=0))
    name_html = input_dir+"/"+naming+"_sunburst.html"
    print(name_html)
    fig.write_html(name_html)
    fig.show()
    return data


# In[399]:


#save_csv = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/unique_final_suspectlist.tsv", sep = "\t")


# In[400]:


#save_csv.to_csv("/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/unique_final_suspectlist.csv")


# # Results analysis

# In[323]:


# 1. list all files


# In[403]:


def list_all_folders(input_dir, list_of_words_to_remove, starting_words_for_all_folders = None):
    #define the input directory
    path = input_dir
    # list all files
    file = os.listdir(path)
    if starting_words_for_all_folders:
        # list all folders that start with "DS" and arent mzML files
        folders = [x for x in os.listdir(path) if x.startswith(starting_words_for_all_folders)]
    else:
        # list all folders that start with "DS" and arent mzML files
        folders = [x for x in os.listdir(path)]
    
    for i in list_of_words_to_remove:
        folders = [x for x in folders if not i in x]
        
    return folders


# In[404]:


input_dir = "./SmarinoiRun1"
list_of_words_to_remove = [".mzML"]
starting_words_for_all_folders = "DS"


# In[405]:


folders2 = list_all_folders(input_dir, list_of_words_to_remove, starting_words_for_all_folders)
folders2


# In[ ]:





# In[ ]:


# 2. Merge all files from MAW and draw a SunBurst


# In[ ]:


# merge all features from all files into one
all_msi = []
for i in folders2:
    entry = path+ "/" +i
    MSI_level_imp = entry+ "/mergedResults-with-one-Candidates.csv"
    if os.path.exists(MSI_level_imp):
        all_msi.append(MSI_level_imp)
        
all_msi_df = pd.concat(map(pd.read_csv, all_msi), ignore_index=True)
all_msi_df
# remove unimportant columns

all_msi_df = all_msi_df[['id_X_x',
       'premz', 'rtmed', 'rtmean', 'int', 'col_eng', 'pol', 'ms2Peaks',
       'ms1Peaks', 'sirius_result_dir', 'id_X_y', 'rtmin', 'rtmax',
       'source_file', 'mbank_results_csv', 'gnps_results_csv',
       'hmdb_results_csv', 'Formula', 'SMILES', 'PubChemID', 'IUPAC',
       'synonyms', 'AnnotationSources', 'AnnotationCount', 'MSILevel', 'MCSS',
       'superclass', 'class', 'subclass', 'ClassificationSource']]
# remove all duplicated features
all_msi_df.drop_duplicates(subset=['premz','rtmed'])
# save the all features from all files into all_features_candidates.csv
all_msi_df.to_csv(path+"/all_features_candidates.csv")


# In[ ]:


# MSI level numbers


# In[565]:


all_msi_df1 = pd.read_csv(path+"/all_features_candidates.csv")
# no. of MSI 2 level
len(all_msi_df1[all_msi_df1["MSILevel"] == 2.0])
# no. of MSI 3 level
len(all_msi_df1[all_msi_df1["MSILevel"] == 3.0])
# other levels
MSI_LEVEL4and5 = all_msi_df1[all_msi_df1["MSILevel"] != 3.0][all_msi_df1["MSILevel"] != 2.0]


# In[ ]:


# add suspect list to analysis


# In[163]:


# list of compounds from SUSPECT LIST
# remove CHOLINE

# 2-O-alpha-L-Mannopyranosyl-D-glyceric acid
# remove Asparagine
# remove Aspartate
# Stearic acid

# remove Succinic acid

# indole
# 2-(2,3-Dihydroxypropoxy)-6-(hydroxymethyl)oxane-3,4,5-triol
# NN-Dimethyl glycine
# to check Butanoic acid
# 2-Methyl-4-carboxy-3,4,5,6-tetrahydropyrimidine
# remove Serine
# Sulfoacetaldehyde
# 3-[3-(3,4-Dihydroxyphenyl)prop-2-enoyloxy]-1,4,5-trihydroxy-cyclohexanecarboxylic acid
# 2,3,4-Trihydroxybutanoic acid
# [1-decanoyloxy-3-[3,4,5-trihydroxy-6-[[3,4,5-trihydroxy-6-(hydroxymethyl)oxan-2-yl]oxymethyl]oxan-2-yl]oxypropan-2-yl] undecanoate  
# Myristic Acid
# 5-(1,2-Dihydroxyethyl)oxolane-2,3,4-trione
# to check D-Xylose ()some industrial interest


# In[ ]:


# suspect list analysis


# In[418]:


suspect = pd.read_csv("/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/unique_final_suspectlist.csv")
suspect = suspect.drop_duplicates(subset= 'SMILES')
# save unique suspect list
suspect.to_csv(
"/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/unique_final_suspectlist.csv")
only_classes = suspect.dropna(subset = 'superclass')
only_classes.to_csv("/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/unique_final_suspectlist_superclasses.csv")


# In[415]:


# plot sunburst again for unique suyspect list ONLY
sunburst(input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList",
         input_csv = "/Users/mahnoorzulfiqar/OneDriveUNI/SuspectList/unique_final_suspectlist_superclasses.csv",
        naming = "suspectlist")


# In[ ]:


# Suspet list comparison


# In[426]:


# read the file with manually added common names
all_msi_df = pd.read_csv(path+"/all_features_candidates_withNames.csv", sep = ";")
# remove SMILE/ compound duplicates
all_msi_df_for_sl = all_msi_df.drop_duplicates(subset= 'SMILES')
# from MAW results and SUSPECT LIST, calculate how many MAW results are same as suspect list:
# y contains the indices of compounds from suspect list which are same as MAW
y = []
# MAW loop
for i, row in all_msi_df_for_sl.iterrows():
    # suspect list loop
    for j, row in suspect.iterrows():
        # if no NAN SMILES
         if not isNaN(all_msi_df_for_sl["SMILES"][i]) and not isNaN(suspect["SMILES"][j]):
                #CALCULATE SMILES
                SKms = [
                    Chem.MolFromSmiles(all_msi_df_for_sl["SMILES"][i]),
                    Chem.MolFromSmiles(suspect["SMILES"][j]),
                ]
                SKfps = [
                    AllChem.GetMorganFingerprintAsBitVect(x, 2, nBits=2048)
                    for x in SKms
                ]
                SKtn = DataStructs.FingerprintSimilarity(SKfps[0], SKfps[1])
                # IF identical SMILES
                if SKtn >= 0.99:
                    # add suspect list index to y
                    y.append(j)
                    print(SKtn)
                    print(suspect["SMILES"][j])
                    # add suspect list to annotation source
                    all_msi_df_for_sl['AnnotationSources'][i] = all_msi_df_for_sl['AnnotationSources'][i] + "|suspectlist"
                    print(suspect["CompoundName"][j])
                    # if couldnt find names for MAW, add names from suspect list
                    if isNaN(all_msi_df_for_sl["Name"][i]):
                        all_msi_df_for_sl["Name"][i] = suspect["CompoundName"][j]


# In[428]:


# x contains indices of compounds from MAW same as 
x = []
for i, rows in all_msi_df_for_sl.iterrows():
    if "suspectlist" in str(all_msi_df_for_sl['AnnotationSources'][i]):
        x.append(i)
        
# intersection
# list of compounds from MAW that also common in suspect list
df2=all_msi_df_for_sl.query('index in @x')
len(df2)


# In[571]:


df2.to_csv(path+"/intersection.csv")


# In[572]:


sunburst(input_dir = path,
        input_csv = path+"/intersection.csv",
        naming = "insection_sl_sm")


# In[429]:


#union


# In[430]:


# MAW without suspect list
MAW_wo_sl = all_msi_df_for_sl.query('index not in @x')
MAW_wo_sl


# In[431]:


# suspect list without MAW
suspect_wo_MAW = suspect.query('index not in @y')
suspect_wo_MAW


# In[432]:


list(suspect_wo_MAW.columns)


# In[450]:


#del suspect_wo_MAW['Unnamed: 0']


# In[434]:


suspect_wo_MAW.columns = ['Name', 
                          'Formula', 
                          'Species', 
                          'SMILES', 
                          'InChI', 
                          'MonoisotopicMass', 
                          'ChEBIid', 
                          'KEGGid', 
                          'PubChemID', 
                          'source_database',
                          'Source', 
                          'nonIsomeric_SMILES_byRDKit', 
                          'IUPAC', 
                          'synonyms',
                          'PubChemPY',
                          'correct_CompoundName',
                          'Molecular mass',
                          'subclass', 
                          'class',
                          'superclass', 
                          'Enzymes', 
                          'InChIKey']


# In[435]:


del MAW_wo_sl["Unnamed: 0"]


# In[436]:


MAW_wo_sl.columns


# In[437]:


suspect_wo_MAW.columns


# In[438]:


all_msi_df_for_sl.append(suspect_wo_MAW, sort=True,ignore_index=True)


# In[439]:


combined_union = all_msi_df_for_sl.append(suspect_wo_MAW, sort=True,ignore_index=True)
combined_union


# In[440]:


combined_union_col = pd.DataFrame(combined_union_col)


# In[441]:


combined_union_col = combined_union.reindex(columns=['Name',
                                                    'SMILES',
                                                    'Formula',
                                                     'Molecular mass',
                                                     'IUPAC',
                                                     'InChI',
                                                     'PubChemID',
                                                    'superclass',
                                                    'class',
                                                    'subclass',
                                                    'synonyms',
                                                     'MSILevel',
                                                     'AnnotationSources',
                                                    'source_file',
                                                    'source_database',
                                                     'Source',
                                                     'Species'
                                                    ])


# In[442]:


combined_union_col.to_csv("./SmarinoiUnion.csv")


# In[443]:


combined_union_col = pd.read_csv("./SmarinoiUnion.csv")


# In[444]:


combined_union_col


# In[445]:


# list of metabolites from other sources than diatom or bacteria
list_met = ["Diflufenican", "L-Homocitrulline", "Succinic Acid", "D-Epiandrosterone", "Captopril", "Succinic semialdehyde",
            "Gabapentin", "Aldimorph", "Adamantane", "Methylsulfocyanat", "p-Octopamine", "Tetradecanamide", "Bursaphelocide B",
            "Paracetamol", "Demethoxycurcumin", "Cyanidin", "8'-Hydroxyrotenone", "Methyl nicotinate", "Rutin", 
            "Bafilomycin D", "Destruxin B", "3'-Demethylnobiletin", "Pomiferin", "Kotanin", "Methyl mycophenolate",
            "Theobromine", "cianidanol", "Puerarin", 'Alpha-Lactose']


# In[446]:


for i in list_met:
    print(i)
    combined_union_col = combined_union_col[combined_union_col.Name != i]


# In[448]:


combined_union_col


# In[449]:


combined_union_col.to_csv("./SmarinoiRun1/SmarinoiUnion_curated.csv")


# In[544]:


input_dir = "./SmarinoiRun1"
input_csv = "./SmarinoiRun1/SmarinoiUnion_curated.csv"
naming = "union"


# In[545]:


sunburst(input_dir, input_csv, naming)


# In[607]:


def chemMN(input_dir, input_csv, naming, name_col):
    df = pd.read_csv(input_csv)
    dbn= []
    for i, row in df.iterrows():
        for j, row in df.iterrows():
            if df['SMILES'][i] != df['SMILES'][j]:
                try:
                    ms = [Chem.MolFromSmiles(df['SMILES'][i]), Chem.MolFromSmiles(df['SMILES'][j])]
                    fps = [AllChem.GetMorganFingerprintAsBitVect(x,2, nBits=1024) for x in ms]
                    tn = DataStructs.FingerprintSimilarity(fps[0],fps[1])
                    dbn.append({
                        'Namei':df[name_col][i],
                        'Namej':df[name_col][j],
                        'i': df['SMILES'][i],
                        'j': df['SMILES'][j],
                        'Tanimoto': tn
                    })
                except:
                    pass
        #print(i)
    db_edge = pd.DataFrame(dbn)
    db_edge.to_csv(input_dir+ "/"+ naming+ "_allVSall.csv")

    dfe = []
    x=0
    for i, row in db_edge.iterrows():        
        if db_edge['Tanimoto'][i] >= 0.85:
            x=x+1
            dfe.append({
                'Start':db_edge['Namei'][i],
                'End':db_edge['Namej'][i],
                'Tanimoto':db_edge['Tanimoto'][i]
            })
    new_df = pd.DataFrame(dfe)
    new_df['Start'] = new_df['Start'].astype(str)
    new_df['End'] = new_df['End'].astype(str)
    new_df['StartAtt']=np.nan
    new_df['EndAtt']=np.nan
    for i, row in new_df.iterrows():
        for j, row in df.iterrows():
            if new_df['Start'][i]==df[name_col][j]:
                new_df.loc[i, 'StartAtt'] = df['superclass'][j]
    for i, row in new_df.iterrows():
        for j, row in df.iterrows():
            if new_df['End'][i]==df[name_col][j]:
                new_df.loc[i, 'EndAtt'] = df['superclass'][j]

    new_df['sorted_names'] = new_df.apply(lambda row: '-'.join(sorted([row['Start'], row['End']])), axis=1)
    new_df = new_df.drop_duplicates(subset=["sorted_names"], keep="last")
    new_df.to_csv(input_dir + "/" + naming + "_chemMN_Cytoscape.tsv", sep='\t')
    return new_df


# In[599]:


chemMN(
       input_dir = "/Users/mahnoorzulfiqar/OneDriveUNI/MAW-Diatom/SmarinoiRun1",
        input_csv = "/Users/mahnoorzulfiqar/OneDriveUNI/MAW-Diatom/SmarinoiRun1/SmarinoiUnion_curated.csv",
       naming = "union",
      name_col = "Name")


# In[ ]:


chemMN(
       input_dir = "./SmarinoiRun1",
        input_csv = "./SmarinoiRun1/SmarinoiUnion_curated.csv",
       naming = "union",
      name_col = "Name")


# In[ ]:


chemMN(
       input_dir = "./suspectlist",
        input_csv = "./suspectlist/unique_final_suspectlist.csv",
       naming = "union",
      name_col = "Name")

