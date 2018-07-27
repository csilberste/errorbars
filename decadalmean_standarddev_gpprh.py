# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 13:44:24 2018

@author: Ceci Silberstein
"""

import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import output_utils as out
from matplotlib.backends.backend_pdf import PdfPages

'''
A potential issue with this function is that if the length of the original array
is not divisible by the grouping interval, the sum of the last remaining elements
of the original array after the for loop will be scaled by a factor of
(grouping)/(length of the remaining array). This means if the array being grouped
had, say, 1310 elements and was being grouped by 12s, the last 2 elements would be
grouped into the resulting array as 6 times their value
MAYBE THIS WILL SOMEDAY BE RELEVANT ONCE AGAIN

def array_grouping(origarray, grouping):
    groupedArray = []
    groupRange = len(origarray)/grouping
    for i in range(0, groupRange):
        newSum = np.sum(origarray[i*grouping:(i+1)*grouping])
        groupedArray.append(newSum)
    if groupRange*grouping != len(origarray):
        print 'GOT HERE!'
        lastSum = np.sum((origarray[groupRange*grouping:])*grouping//len(origarray[groupRange*grouping:]))
        groupedArray.append(lastSum)
    return groupedArray
'''
    

'''
'''
def separate_means(origarray, meanInterval):
    meanARR = []
    meanRange = len(origarray)//meanInterval
    #print meanRange
    #print 'meanrange'
    for i in range(0, meanRange):
       # print (i*meanRange)
       # print ((1+i)*meanRange)
       # print  np.mean(origarray[(i*meanRange):((1+i)*meanRange)])
        localMean = np.mean(origarray[(i*meanInterval):((1+i)*meanInterval)])
        meanARR.append(localMean)
        
    if len(origarray) - meanInterval*meanRange > 7:
        lastMean = np.mean(origarray[(meanInterval*meanRange):])
        meanARR.append(lastMean)
  #  print meanARR
    return meanARR


'''
 
'''
def standard_deviation_multipiece(origarray, meanarr):
    interval = len(origarray)/len(meanarr)
    standarddevARR = []
    '''
    REPLACE/COMPARE with build in numpy function
    '''
    for i in range(0, len(meanarr)):
        sumOfDistSq = 0
        for j in origarray[i*interval:(i+1)*interval]:
            sumOfDistSq +=(j - meanarr[i])**2
        meanOfThat = sumOfDistSq/interval
        sqrtOfMean = meanOfThat**(.5)
        standarddevARR.append(sqrtOfMean)
    #print standarddevARR
    return standarddevARR

'''
This sorry old function takes a list with a length not divisible by 10
and cuts off anything past point where it can be grouped by 10
'''
def tens_truncate(sadness):
    sad_list = []
    sad_list = sadness
    if len(sad_list)%10 == 0:
        #list is NOT SAD
        return sad_list
    else:
        return sad_list[:(len(sad_list)//10)*10]

     
'''
PROCESS
'''
sitesARR = ['1', '2', '3', '4', '5']
'''
Process for gathering arrays of historic, projected climate data
'''
histclimateARR = []
projclimateARR = []
histenvARR = []
projenvARR = []
histyearlyARR = []
projyearlyARR = []
environmental_variables = ['GPP', 'RH']
for env in environmental_variables:    
    for num in sitesARR:
        historicclimate = nc.Dataset('C:/Users/Public/Anaconda/dalton_transect_sites/dalton_transect_sites/site_'+num+'_1x1/historic-climate.nc')
        hclimate = historicclimate.variables['precip'][:,0,0]
        hcgroupARR = out.sum_monthly_flux_to_yearly(hclimate)
        hcmeanARR = separate_means(hcgroupARR, 10)
        hcstanddevs = standard_deviation_multipiece(hcgroupARR, hcmeanARR)
        histclimateARR.append([hcmeanARR, hcstanddevs])
        
        projectedclimate = nc.Dataset('C:/Users/Public/Anaconda/dalton_transect_sites/dalton_transect_sites/site_'+num+'_1x1/projected-climate.nc')
        pclimate = projectedclimate.variables['precip'][:,0,0]
        pcgroupARR = out.sum_monthly_flux_to_yearly(pclimate)
        pcmeanARR = separate_means(pcgroupARR, 10)
        pcstanddevs = standard_deviation_multipiece(pcgroupARR, pcmeanARR)
        projclimateARR.append([pcmeanARR, pcstanddevs])
        
        historicenv = nc.Dataset('C:/Users/Public/Anaconda/Runs_ALLPFT/DHS'+num+'_CMT04/'+env+'_monthly_tr.nc')
        henv = historicenv.variables[env][:,0,0]
        henvgroupARR = out.sum_monthly_flux_to_yearly(henv)
        histyearlyARR.append(henvgroupARR)
        henvmeanARR = separate_means(henvgroupARR, 10)
        henvstanddevs = standard_deviation_multipiece(henvgroupARR, henvmeanARR)
        histenvARR.append([henvmeanARR, henvstanddevs])
        
        projectedenv= nc.Dataset('C:/Users/Public/Anaconda/Runs_ALLPFT/DHS'+num+'_CMT04/'+env+'_monthly_sc.nc')
        penvt = projectedenv.variables[env][:,0,0]
        penvgroupARR = out.sum_monthly_flux_to_yearly(penvt)
        projyearlyARR.append(penvgroupARR)
        penvmeanARR = separate_means(penvgroupARR, 10)
        penvstanddevs = standard_deviation_multipiece(penvgroupARR, penvmeanARR)
        projenvARR.append([penvmeanARR, penvstanddevs])
        
for num in range(0, len(sitesARR)):
    henvgroupARR = np.subtract(histyearlyARR[num], histyearlyARR[num+len(sitesARR)])
    henvmeanARR = separate_means(henvgroupARR, 10)
    henvstanddevs = standard_deviation_multipiece(henvgroupARR, henvmeanARR)
    histenvARR.append([henvmeanARR, henvstanddevs])
    
    penvgroupARR = np.subtract(projyearlyARR[num], projyearlyARR[num+len(sitesARR)])
    penvmeanARR = separate_means(penvgroupARR, 10)
    penvstanddevs = standard_deviation_multipiece(penvgroupARR, penvmeanARR)
    projenvARR.append([penvmeanARR, penvstanddevs])
        

'''
First attempt at error bar plot
'''
for env in environmental_variables:
    pdf = PdfPages('C:/Users/Public/Anaconda/dalton_transect_sites/'+env+'ERRORBARS.pdf')
    
    for i in range(0, len(histenvARR[0][0])):
        fig, ax = plt.subplots()
        for j in range(0, len(sitesARR)):
            if env == 'RH':
                k = j+len(sitesARR)
            else:
                k = j
            ax.errorbar(histclimateARR[k][0][i], histenvARR[k][0][i], histenvARR[k][1][i], histclimateARR[k][1][i], label=sitesARR[j])
        ax.set_xlabel('Annual Precipitation')
        if env == 'GPP':
            ax.set_ylabel('Gross Primary Productivity')
        else:
            ax.set_ylabel('Heterotrophic Respiration')
        ax.set_title(str(1900+ 10*i))
        if env == 'GPP':
            ax.set_xlim(100, 900)
            ax.set_ylim(350, 1000)
        else:
            ax.set_xlim(100, 900)
            ax.set_ylim(-300, 400)
        ax.legend()
        pdf.savefig()
    
    for i in range(0, len(projenvARR[0][0])):
        fig, ax = plt.subplots()
        for j in range(0, len(sitesARR)):
            if env == 'RH':
                k = j+len(sitesARR)
            else:
                k = j
            ax.errorbar(projclimateARR[k][0][i], projenvARR[k][0][i], projenvARR[k][1][i], projclimateARR[k][1][i], label=sitesARR[j])
        ax.set_xlabel('Annual Precipitation')
        if env == 'GPP':
            ax.set_ylabel('Gross Primary Productivity')
        else:
            ax.set_ylabel('Heterotrophic Respiration')
        ax.set_title(str(2010+ 10*i))
        if env == 'GPP':
            ax.set_xlim(100, 900)
            ax.set_ylim(350, 1000)
        else:
            ax.set_xlim(100, 900)
            ax.set_ylim(-300, 400)
        ax.legend()
        pdf.savefig()
    pdf.close()

pdf = PdfPages('C:/Users/Public/Anaconda/dalton_transect_sites/NEPERRORBARS.pdf')
    
for i in range(0, len(histenvARR[0][0])):
    fig, ax = plt.subplots()
    for j in range(0, len(sitesARR)):
        k = j + 2*len(sitesARR)
        ax.errorbar(histclimateARR[j][0][i], histenvARR[k][0][i], histenvARR[k][1][i], histclimateARR[j][1][i], label=sitesARR[j])
    ax.set_xlabel('Annual Precipitation')
    ax.set_ylabel('Net Ecosystem Production')
    ax.set_title(str(1900+ 10*i))
    ax.set_xlim(100, 900)
    ax.set_ylim(0, 1000)
    ax.legend()
    pdf.savefig()
for i in range(0, len(projenvARR[0][0])):
    fig, ax = plt.subplots()
    for j in range(0, len(sitesARR)):
        k = j + 2*len(sitesARR)
        ax.errorbar(projclimateARR[j][0][i], projenvARR[k][0][i], projenvARR[k][1][i], projclimateARR[j][1][i], label=sitesARR[j])
    ax.set_xlabel('Annual Precipitation')
    ax.set_ylabel('Net Ecosystem Production')
    ax.set_title(str(2010+ 10*i))
    ax.set_xlim(100, 900)
    ax.set_ylim(0, 1000)
    ax.legend()
    pdf.savefig()
pdf.close()
    