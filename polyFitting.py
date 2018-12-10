import sys
import os
import csv
import numpy as np
import pandas as pd
import argparse

def search(lib,row):
    if row['INCHI'] in lib:
        return lib[row['INCHI']]
    return -1

def main():
    # parse the argument
    parser = argparse.ArgumentParser(description='run polyFitting')
    parser.add_argument('input', help='input')
    parser.add_argument('cosineScore', help='cosineScore')
    parser.add_argument('error_threshold', help='how much the\
                        error sets to be to define outliers')
    args = parser.parse_args()
    inputF= args.input
    cosineScore = float(args.cosineScore)
    error = float(args.error_threshold)

    df = pd.read_csv(inputF,sep='\t')
    # filter out by cosinScore
    df = df[df.MQScore > cosineScore]
    # search for ki average
    lib_df = pd.read_csv('nonpolar.csv')
    lib_df = lib_df[lib_df.polarity.str.contains('non-polar')]
    lib = pd.Series(lib_df.ki_nonpolar_average.values,index=lib_df.INCHI.values).to_dict()
    df['ki_average'] = df.apply(lambda row:search(lib,row),axis = 1)
    df = df[df.ki_average>0]
    
    #clean the data for polynomial fitting:
    df = df[df['ki_average']>500]
    df = df[df.RT_Query<800]
    #simply find the polynomial fitting
    p_a = np.polyfit(df.RT_Query,df.ki_average,2)
    # we fit it twice to have more robust results:
    df =df[abs(df.ki_average-np.polyval(p_a,df.RT_Query))/np.polyval(p_a,df.RT_Query)<error]
    p_b = np.polyfit(df.RT_Query,df.ki_average,2)
    outputF = "polyFittingParameters.txt"
    np.savetxt(outputF, p_b, delimiter=',')

if __name__=="__main__":
    main()
