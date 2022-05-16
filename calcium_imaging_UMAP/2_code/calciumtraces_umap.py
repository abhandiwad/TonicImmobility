# -*- coding: utf-8 -*-
"""
Created on Mon Feb  7 09:32:49 2022

This script:
    
    1. Takes in normalized calcium imaging traces
    2. Applies the continuous wavelet transform to each trace
    3. Outputs a UMAP figure that clusters cells based on wavelet-transformed calcium signal dynamics


@author: Bhandiwad
"""

import os
import pywt
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import pandas as pd
import umap

def read_calcium_traces():
    
    catraces = pd.read_table('Fig7E_calciumtraces.csv',sep=',')
    catraces = catraces.to_numpy()
    
    return catraces

def plot_wavelet(time, signal, scales, 
                 waveletname = 'cmor', 
                 title = 'Wavelet Transform (Power Spectrum) of signal', 
                 ylabel = 'Frequency (Hz)', 
                 xlabel = 'Time'):
    
    dt = time[1] - time[0]
    [coefficients, frequencies] = pywt.cwt(signal, scales, waveletname, dt)
    power = (abs(coefficients)) ** 2
    period = 1. / frequencies
    levels = [0.0625, 0.125, 0.25, 0.5, 1, 2, 4, 8]
    contourlevels = np.log2(levels)
    
    fig, ax = plt.subplots(figsize=(15, 10))
    im = ax.contourf(time, np.log2(period), np.log2(power), contourlevels, extend='both')
    
    ax.set_title(title, fontsize=20)
    ax.set_ylabel(ylabel, fontsize=18)
    ax.set_xlabel(xlabel, fontsize=18)
    
    yticks = 2**np.arange(np.ceil(np.log2(period.min())), np.ceil(np.log2(period.max())))
    ax.set_yticks(np.log2(yticks))
    ax.set_yticklabels(yticks)
    ax.invert_yaxis()
    ylim = ax.get_ylim()
    ax.set_ylim(ylim[0], -1)
    
    
    
    cbar_ax = fig.add_axes([0.95, 0.5, 0.03, 0.25])
    fig.colorbar(im, cax=cbar_ax, orientation="vertical")
    sns.set_palette('Reds')
    plt.show()
    
    
def wavelet_power(time, signal, scales, waveletname = 'cmor'):
    dt = time[1] - time[0]
    [coefficients, frequencies] = pywt.cwt(signal, scales, waveletname, dt)
    power = (abs(coefficients)) ** 2
    
    #plot_wavelet(time, signal, scales)
    
    return power

def wavelet_transform_all(catraces):
    
    time = catraces[:,0]
    scales = np.arange(1,128) # Make sure this value less less than the Nyquist frequency

    wavetransformd = np.empty((np.shape(catraces)[1],len(scales)*np.shape(catraces)[0]))

    for i in range(1,np.shape(catraces)[1]):
        ca = catraces[:,i]
        power = wavelet_power(time,ca,scales)
        
        wavetransformd[i,:] = np.reshape(power,len(scales)*np.shape(catraces)[0])
        
    return wavetransformd

def generate_UMAP(wavetransformd):
    
    wavelets = wavetransformd[1::,:]

    reducer = umap.UMAP(random_state=42)
    reducer.fit(wavelets)

    embedding = reducer.transform(wavelets)
    assert(np.all(embedding == reducer.embedding_))
    #embedding.shape
    
    return embedding
    
#%%
if __name__ == "__main__": 
    os.chdir('D:/Manuscripts/SEBA paper/Figures_data/Figure7')
    
    catraces = read_calcium_traces()
    
    wavetransformd = wavelet_transform_all(catraces)
            
    embedding = generate_UMAP(wavetransformd)   
    
    
    #%% Plot UMAP figure
    
    plt.subplots(figsize=(15, 10))
    plt.scatter(embedding[:, 0],embedding[:, 1])
    plt.gca().set_aspect('equal', 'datalim')
    plt.show()



















