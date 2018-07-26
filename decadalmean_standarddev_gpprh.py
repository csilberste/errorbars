# -*- coding: utf-8 -*-
"""
Created on Thu Jul 12 13:44:24 2018

@author: Rubus occidentalis
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
'''
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
def separate_means(origarray, meanInterval):
    meanARR = []
    meanRange = len(origarray)//meanInterval
    #print meanRange
    #print 'meanrange'
    for i in range(0, meanInterval):
       # print (i*meanRange)
       # print ((1+i)*meanRange)
       # print  np.mean(origarray[(i*meanRange):((1+i)*meanRange)])
        localMean = np.mean(origarray[(i*meanRange):((1+i)*meanRange)])
        meanARR.append(localMean)
    if meanInterval*meanRange != len(origarray):
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
PROCESS
'''
sitesARR = ['1', '2', '3', '4', '5']
'''
Process for gathering arrays of historic, projected climate data
'''
histclimateARR = []
projclimateARR = []
histgppARR = []
projgppARR = []
environmental_variables = ['GPP', 'RH']
for env in environmental_variables:    
    for num in sitesARR:
        historicclimate = nc.Dataset('C:/Users/Public/Anaconda/dalton_transect_sites/dalton_transect_sites/site_'+num+'_1x1/historic-climate.nc')
        hclimate = historicclimate.variables['precip'][:,0,0]
        hclimateARR = ma.getdata(hclimate)
        hcgroupARR = array_grouping(hclimateARR, 12)
        hcmeanARR = separate_means(hcgroupARR, 10)
        hcstanddevs = standard_deviation_multipiece(hcgroupARR, hcmeanARR)
        histclimateARR.append([hcmeanARR, hcstanddevs])
        
        projectedclimate = nc.Dataset('C:/Users/Public/Anaconda/dalton_transect_sites/dalton_transect_sites/site_'+num+'_1x1/projected-climate.nc')
        pclimate = projectedclimate.variables['precip'][:,0,0]
        pclimateARR = ma.getdata(pclimate)
        pcgroupARR = array_grouping(pclimateARR, 12)
        pcmeanARR = separate_means(pcgroupARR, 10)
        pcstanddevs = standard_deviation_multipiece(pcgroupARR, pcmeanARR)
        projclimateARR.append([pcmeanARR, pcstanddevs])
        
        historicgpp = nc.Dataset('C:/Users/Public/Anaconda/Runs_ALLPFT/DHS'+num+'_CMT04/'+env+'_monthly_tr.nc')
        hgpp = historicgpp.variables[env][:,0,0]
        hgppARR = ma.getdata(hgpp)
        hgppgroupARR = array_grouping(hgppARR, 12)
        hgppmeanARR = separate_means(hgppgroupARR, 10)
        hgppstanddevs = standard_deviation_multipiece(hgppgroupARR, hgppmeanARR)
        histgppARR.append([hgppmeanARR, hgppstanddevs])
        
        projectedgpp= nc.Dataset('C:/Users/Public/Anaconda/Runs_ALLPFT/DHS'+num+'_CMT04/'+env+'_monthly_sc.nc')
        pgpp = projectedgpp.variables[env][:,0,0]
        pgppARR = ma.getdata(pgpp)
        pgppgroupARR = array_grouping(pgppARR, 12)
        pgppmeanARR = separate_means(pgppgroupARR, 10)
        pgppstanddevs = standard_deviation_multipiece(pgppgroupARR, pgppmeanARR)
        projgppARR.append([pgppmeanARR, pgppstanddevs])

'''
First attempt at error bar plot

for env in environmental_variables:
    pdf = PdfPages('C:/Users/Public/Anaconda/dalton_transect_sites/'+env+'ERRORBARS.pdf')
    for i in range(0, len(histgppARR[0][0])):
        fig, ax = plt.subplots()
        for j in range(0, len(sitesARR)):
            if env == 'RH':
                k = j+len(sitesARR)
            else:
                k = j
            ax.errorbar(histclimateARR[k][0][i], histgppARR[k][0][i], projgppARR[k][1][i], projclimateARR[k][1][i], label=sitesARR[j])
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
    for i in range(0, len(projgppARR[0][0])):
        fig, ax = plt.subplots()
        for j in range(0, len(sitesARR)):
            if env == 'RH':
                k = j+len(sitesARR)
            else:
                k = j
            ax.errorbar(projclimateARR[k][0][i], projgppARR[k][0][i], projgppARR[k][1][i], projclimateARR[k][1][i], label=sitesARR[j])
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
'''
print len(projgppARR[0][0])