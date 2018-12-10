# -*- coding: utf-8 -*-
"""Kovat Index Filter
This file is used to filter out the compounds with bad retention Queries.

Example
-------
    $ python mapping.py input.csv mode(p/m) polyfitParameter/carbonMarker.csv cosineScore(float) errorTolerance(float)


Notes
-----
Run it with kovatLib.csv in the same directory.
"""

import sys
import os
import csv 
import numpy as np
import pandas as pd
import argparse



# load df from input data file
def loadDf(csv,cosine,prediction,mode):
    df = pd.read_csv(csv,sep = '\t')
    new_df = pd.DataFrame({'#Scan#':df['#Scan#'],	
                           'CAS': df['CAS_Number'],'Name':df['Compound_Name'],
                           'Cosine':df['MQScore'], 'INCHI':df['INCHI'],
                           'ki_estimate':np.nan,'ki_average':np.nan,
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
        #drop unidentified Items
        if new_df['CAS'][i]  == 'N/A'and new_df['Name'][i] == 'N/A':
            toDrop.append(i)
    new_df = new_df.drop(i for i in toDrop)
    new_df = new_df.reset_index(drop=True)

    # load database:
    # library is sorted by CAS
    lib_df = pd.read_csv('nonpolar.csv')
    lib_df = lib_df[lib_df.polarity.str.contains('non-polar')]
    lib = pd.Series(lib_df.ki_nonpolar_average.values,index=lib_df.INCHI.values).to_dict()


    #fill in the kovat index from the library search
    for i in range(len(new_df)):
        try:
            new_df['ki_average'][i] = lib[new_df['INCHI'][i]]
            if mode == 'm':
                new_df['ki_estimate'][i] =kovatIndex(float(new_df['RT'][i]), prediction)
            else:
                new_df['ki_estimate'][i]= np.polyval(prediction,float(new_df['RT'][i]))
            new_df['Error'][i] = abs(new_df['ki_estimate'][i] - new_df['ki_average'][i])/new_df['ki_average'][i]
        except:
            continue       
    return new_df

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
            return 0
        elif (rt > markerDic.RT_Query[i] and rt < markerDic.RT_Query[i+1]) \
             or (rt == markerDic.RT_Query[i] or rt ==  markerDic.RT_Query[i+1]):
            N,n,tr_N,tr_n = markerDic['Compound_Name'][i+1],markerDic['Compound_Name'][i], \
                           markerDic.RT_Query[i+1],markerDic.RT_Query[i]
            ki_estimate = 100.0*(n+(N-n)*(rt-tr_n)/(tr_N - tr_n))
            return ki_estimate
    return 0.0
    
    
def main():
    # parse the argument
    parser = argparse.ArgumentParser(description='run mapping')
    parser.add_argument('input', help='input')
    parser.add_argument('mode',help='polynomial fitting or carbon marker enter (p/m)')
    parser.add_argument('additionalFile', help='additional file')
    parser.add_argument('cosineScore', help='cosineScore')
    parser.add_argument('errorTolerance', help='errorTolerance')
    
    args = parser.parse_args()
    inputF= args.input
    mode = args.mode
    additionalFile = args.additionalFile
    cosineScore = float(args.cosineScore)
    errorTolerance = float(args.errorTolerance)

    # load markers
    if mode == 'p':
        prediction = np.loadtxt(additionalFile)
    else:
        prediction = loadMarkers(additionalFile)

    #LoadDataframe and trim it to be the one we need
    df = loadDf(inputF,cosineScore,prediction,mode)
    df.to_csv('nonfiltered.tsv', sep='\t') 

    #filter out the data with the error bigger than the tolerance
    df = df[df.Error < errorTolerance]
    df.to_csv('output.tsv', sep='\t') 
        

if __name__=="__main__":
    main()    
