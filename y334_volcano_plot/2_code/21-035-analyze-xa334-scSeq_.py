# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 12:48:11 2021

@author: Ashwin

This script: 
    
    1. Takes in csv of normalized read counts
    2. Compares xa334 samples to non-vPPN samples
    3. Plots a half-volcano plot with significant genes labeled (Update: this is Figure 4H in paper, 12/27/21 AB)
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os,sys
import scipy.stats as stats
import statsmodels.stats.multitest as smm
import matplotlib.ticker as ticker

thisScript = os.path.basename( sys.argv[0]).split('_')[0]

def set_plotting_params(): 
    
    '''Sets plotting aesthetics'''
    global graphLabelSize
    graphLabelSize = 10
    
    plt.rcParams['font.sans-serif'] = "Arial"
    plt.rcParams['font.family'] = "sans-serif"
    sns.set_style('white', {'xtick.bottom': True,'ytick.left': True })
    sns.set_style( {'axes.spines.right': False, 'axes.spines.top': False } )
    plt.rc('xtick', labelsize=graphLabelSize)
    plt.rc('ytick', labelsize=graphLabelSize)
    

    

def read_data_file(inFile): 
    '''Reads in csv containing normalized counts'''

    df = pd.read_csv(  inFile , index_col=0 ,sep=' ')
    df = df.set_index('gene_name')
    
    return df

def remove_zeros(df):  
    
    '''remove points with no reads in any column'''
    df['m']=df.apply(lambda x: max(x), axis=1)
    df = df[ df['m'] > 0 ] 
    df = df.drop('m', axis=1)
    
    return df

def median_counts(df):
    '''get medians of remaining columns'''
    df['ctrl'] = df[ ['elavA', 'elavB', 'elavC', 'elavD', 'elavE', 'elavF', 'elavG', 'elavH', 'elavI', 'x330A', 'x330B', 'x330C', 'x104A', 'x104B', 'x104C', 'x240A', 'x240B', 'x240C', 'x240D','x108A','x108B','x108C'] ].median(axis=1)
    df['x334'] = df[ ['x334A', 'x334B', 'x334C', 'x334D','x334E', 'x334F'] ].median(axis=1)
    
    return df

def calculate_log2FC(df):
    
    # filter for genes highly expressed in xa334
    df = df[ df['x334'] > 5000 ]   

    # Remove zero reads and replace with 1 to make log transform viable
    removeZero = lambda x: max(x,1) # for stats, eliminate zero reads - make 1 read.
    df = df.applymap( removeZero )
    df['ctrl'] = df['ctrl'].apply(lambda x: max(x,1) ) # remove zeros for the ratio
    df['x334'] = df['x334'].apply(lambda x: max(x,1) )

    df['ratio'] =  df['x334'] / df['ctrl']
    df['log2ratio'] = np.log2( df['ratio'] )
    
    return df

def half_volcano_MWU(df): #Perform nonparametric Mann-Whitney U test to test for significance
    
    # filter for genes enriched in xa334
    df = df[ df['ratio'] > 1]

    # calculate p values
    df['pval'] =  df.apply( lambda x: stats.mannwhitneyu( df.loc[:,('x334A', 'x334B', 'x334C', 'x334D', 'x334E', 'x334F')],
                                                           df.loc[:,('elavA','elavB','elavC','elavD','elavE','elavF','elavG',
                                                                     'elavH','elavI','x330A','x330B','x330C',
                                                                     'x104A','x104B','x104C',
                                                                     'x240A','x240B','x240C','x240D',
                                                                     'x108A','x108B','x108C')],axis=1)[1])


    # Perform Benjamini-Hochberg correction for FDR
    bh_corr = smm.multipletests(df['pval'],alpha=0.1,method='fdr_bh')
    df['pval_corr'] = bh_corr[1]
    
    # Log transform for volcano plot
    df['log10pval'] = -1*np.log10( df['pval_corr'] )
    
    return df


def plot_half_volcano(df): # half volcano plot for genes enriched in xa334
    
    plt.figure(figsize=(1.7,2.8)) 
    ax = sns.scatterplot(x='log2ratio',y='log10pval', data=df,s=10,alpha=0.5)
    
    # Draw horizontal line marking statistical significance
    padj = -1*np.log10( 0.1  )
    ax.plot( [0,15], [padj,padj], color='red', zorder=0, ls=':' ,lw=0.5)
    
    
    # Label significant genes
    anlist = list( df[ (df['log2ratio'] > 10) & (df['pval_corr'] < 0.1) ].index )
    for an in anlist:
        plt.text( df.loc[an,'log2ratio'], df.loc[an,'log10pval'],an, ha='right',size=graphLabelSize)

    
    ax.set_ylabel('log10(p-value)',size=graphLabelSize, fontweight='bold')
    ax.set_xlabel('log2FC',size=graphLabelSize, fontweight='bold')
    ax.set(ylim=(0,3) )

    ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
    ax.yaxis.set_major_locator(ticker.MultipleLocator(1))
    plt.tight_layout()

#%%

if __name__ == "__main__": 
    
    set_plotting_params()
    
    
    df = read_data_file('C:\\Users\\optim\\Desktop\\ArrestManuscript\\OngoingFigureData\\RNAseq\\normalized_counts_neurons-smartseq.csv')
    
    df = remove_zeros(df)
    df = median_counts(df)
    df = calculate_log2FC(df)
    
    
    df = df[ df['ratio'] > 1] # Only select for enriched genes for half-volcano
    
    statistical_test = stats.mannwhitneyu( df.loc[:,('x334A', 'x334B', 'x334C', 'x334D', 'x334E', 'x334F')] ,
                                                            df.loc[:,('elavA','elavB','elavC','elavD','elavE','elavF','elavG',
                                                                      'elavH','elavI','x330A','x330B','x330C',
                                                                      'x104A','x104B','x104C',
                                                                      'x240A','x240B','x240C','x240D',
                                                                      'x108A','x108B','x108C')],axis=(1))[1]
    
    df = half_volcano_MWU(df)
    
    plot_half_volcano(df)



