import sys
import os
import csv 
import numpy as np
import pandas as pd
import argparse


# load df from input data file
def loadDf(csv,cosine,markerDic,lib):
    df = pd.read_csv(csv,sep = '\t')
    new_df = pd.DataFrame({'CAS': df['CAS_Number'],'Name':df['Compound_Name'],
                           'Cosine':df['MQScore'], 'INCHI':df['INCHI'],
                            'ki_real':np.nan,'ki_estimate':np.nan,'ki_average':np.nan,
                           'TIC': df['TIC_Query'],'RT':df['RT_Query'],'Error':np.nan})
    
    # filter out by cosine score first
    new_df = new_df[new_df.Cosine > cosine]
    new_df = new_df.reset_index(drop=True)
    new_df = new_df.sort_values(['INCHI','TIC','Cosine'], ascending=[True,False,False]) 
    # restructure the dataframe by deleting repeated elements
    preINCHI = 'N/A'
    toDrop =[]
    for i in range(len(new_df)):
        if  new_df['INCHI'][i] == preINCHI:
            toDrop.append(i)
            continue
        preINCHI = new_df['INCHI'][i]
        #drop unIdentified Item
        if new_df['CAS'][i]  == 'N/A'and new_df['Name'][i] == 'N/A':
            toDrop.append(i)
    new_df = new_df.drop(i for i in toDrop)
    new_df = new_df.reset_index(drop=True)

    #fill in the kovat index from the library search
    for i in range(len(new_df)):
        print(float(i)/float(len(new_df)))
        kiList = libSearch(lib,new_df['CAS'][i],new_df['Name'][i])
        if len(kiList) != 0:           
            new_df['ki_real'][i] = ";".join(kiList)
            new_df['ki_estimate'][i] =kovatIndex(float(new_df['RT'][i]), markerDic)
            new_df['ki_average'][i] = sum(float(i) for i in kiList)/float(len(kiList))
            new_df['Error'][i] = abs(new_df['ki_estimate'][i] - new_df['ki_average'][i])/new_df['ki_average'][i]       
    return new_df


# do a general search with name and CAS numbers
def libSearch (lib,CAS, Name):
    kis = []
    enter = False
    for i in range(len(lib)):
        if 'non-polar' in lib['polarity'][i]:
            # taking too long
            if lib['name'][i] ==  lib['name'][i]:
                if lib['name'][i] == Name:
                    kis.append(str(lib['ki'][i]))
            if lib['CAS numbers'][i] == lib['CAS numbers'][i] and CAS == CAS:
                if lib['CAS numbers'][i] == CAS.replace("-", ""):
                    kis.append(str(lib['ki'][i]))
                    enter =True
                elif enter:
                    break
    return kis

def loadMarkers(marker):
    df = pd.read_csv(marker,sep = ';')
    # compounds name has to be in the format of "name(C#)"
    for i in range(len(df)):
        c_n = (df['Compound_Name'][i].split('(C')[-1]).split(")")[0]
        df.ix[i, 'Compound_Name'] = float(c_n)
        df.ix[i, 'RT_Query'] = float(df['RT_Query'][i])
    df = df.sort_values(['Compound_Name'], ascending=[True])
    return df

def kovatIndex(rt, markerDic):
    for i in range(len(markerDic)):
        if i == len(markerDic)-1:
            print(rt,markerDic.RT_Query[i],markerDic['Compound_Name'][i])
            return 0
        elif (rt > markerDic.RT_Query[i] and rt < markerDic.RT_Query[i+1]) \
             or (rt == markerDic.RT_Query[i] or rt ==  markerDic.RT_Query[i+1]):
            N,n,tr_N,tr_n = markerDic['Compound_Name'][i+1],markerDic['Compound_Name'][i], \
                           markerDic.RT_Query[i+1],markerDic.RT_Query[i]
            print(N,n,tr_N,tr_n,rt)
            ki_estimate = 100.0*(n+(N-n)*(rt-tr_n)/(tr_N - tr_n))
            return ki_estimate
             
    print(rt,markerDic.RT_Query[i],markerDic['Compound_Name'][i])
    return 0.0
    
    
def main():
    # This file is used to filter out the compounds with bad retention Queries
    # Run it with kovatLib.csv in the same directory.
    # python mapping.py input.csv carbonMarker.csv cosineScore(float) errorTolerance(float)
    # parse the argument
    parser = argparse.ArgumentParser(description='run mapping')
    parser.add_argument('input', help='input')
    parser.add_argument('carbonMarker', help='carbonMarker')
    parser.add_argument('cosineScore', help='cosineScore')
    parser.add_argument('errorTolerance', help='errorTolerance')
    args = parser.parse_args()
    inputF= args.input
    marker = args.carbonMarker
    cosineScore = float(args.cosineScore)
    errorTolerance = float(args.errorTolerance)

    # load library 
    # library is sorted by CAS
    lib = pd.read_csv('Kovatlib.csv')

    # load markers
    markerDic = loadMarkers(marker)
    print( markerDic)


    #LoadDataframe and trim it to be the one we need
    df = loadDf(inputF,cosineScore,markerDic,lib)
    df.to_csv('nonfiltered.tsv', sep='\t') 

    #filter out the data with the error bigger than the tolerance
    df = df[df.Error < errorTolerance  ]
    df.to_csv('output.tsv', sep='\t') 


    


        

if __name__=="__main__":
    main()    
