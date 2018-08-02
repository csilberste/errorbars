# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 13:19:00 2018

@author: Rubus occidentalis
"""

import os
import numpy as np
import numpy.ma as ma
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import netCDF4 as nc
import output_utils as out
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy import stats
#from sklearn import datasets, linear_model
#from sklearn.metrics import mean_squared_error, r2_score

'''
This function takes divides a list of values into chunks the length of the
meanInterval and returns a list of the means of those chunks

@param list origarray
@param int meanInterval
'''
def separate_means(origarray, meanInterval):
    meanARR = []
    meanRange = len(origarray)//meanInterval
    for i in range(0, meanRange):
        localMean = np.mean(origarray[(i*meanInterval):((1+i)*meanInterval)])
        meanARR.append(localMean)        
    if len(origarray) - meanInterval*meanRange > 7:
        lastMean = np.mean(origarray[(meanInterval*meanRange):])
        meanARR.append(lastMean)
    return meanARR

'''
This function takes a list of values and a list of the means of chunks of those
values and calculates the standard deviation by chunk
**note: In its current form this is done manually, though I intend to replace the
process with numpy.std
****response to note: since "separate_means" sorts for relevant vs irrelevant
uneven chunks at the end of the list (after the portion divisible into 10s), it
seems the built-in standard deviation function will not handle means as I intend
to. This is likely resolvable by running that check for the uneven chunks in the
same single function. For now, however, I intend to focus on functionality before
returning to make it sleeker

@param list origarray
@param list meanarr
'''
def standard_deviation_multipiece(origarray, meanarr):
    interval = len(origarray)/len(meanarr)
    standarddevARR = []
    for i in range(0, len(meanarr)):
        sumOfDistSq = 0
        for j in origarray[i*interval:(i+1)*interval]:
            sumOfDistSq +=(j - meanarr[i])**2
        meanOfThat = sumOfDistSq/interval
        sqrtOfMean = meanOfThat**(.5)
        standarddevARR.append(sqrtOfMean)
    return standarddevARR
'''
This function takes a filename (which should be a netCDF input or output of
dvm-dos-tem with monthly data) and the associated variable (env) and returns the means
of chunks of the list separated by the length of mean_interval.

@param Path filename
@param String env
@param int mean_interval
'''
def file_to_mean_standdev(filename, env, mean_interval):
        whole_file = nc.Dataset(filename)
        maskedarray_fromfile = whole_file.variables[env][:,0,0]
        if env == 'tair':
            yearly_data = out.average_monthly_pool_to_yearly(maskedarray_fromfile)
            #yearly_data = out.sum_monthly_flux_to_yearly(maskedarray_fromfile)
            #print 'got here!'
        else:
            yearly_data = out.sum_monthly_flux_to_yearly(maskedarray_fromfile)
        mean_list = separate_means(yearly_data, mean_interval)
        standdev_list = standard_deviation_multipiece(yearly_data, mean_list)
        return[mean_list, standdev_list]
        
def file_to_mean_yr(filename, env, mean_interval):
        whole_file = nc.Dataset(filename)
        maskedarray_fromfile = whole_file.variables[env][:,0,0]
        unit = whole_file.variables[env].units
        if env == 'tair':
            yearly_data = out.average_monthly_pool_to_yearly(maskedarray_fromfile)
            #yearly_data = out.sum_monthly_flux_to_yearly(maskedarray_fromfile)
            #print 'got here!'
        else:
            yearly_data = out.sum_monthly_flux_to_yearly(maskedarray_fromfile)
        mean_list = separate_means(yearly_data, mean_interval)
        #standdev_list = standard_deviation_multipiece(yearly_data, mean_list)
        return[mean_list, yearly_data, unit]

sitesARR = ['1', '2', '3', '4', '5']
climate_variables = ['tair', 'precip']
environmental_variables = ['NPP', 'RH']
community_type = ['CMT04', 'CMT05', 'CMT06']
climate_dir = 'C:/Users/Public/Anaconda/dalton_transect_sites/dalton_transect_sites/'
env_dir = 'C:/Users/Public/Anaconda/Runs_ALLPFT/'
timeperiod = [['projected', 'sc']] #['historic', 'tr'], 
length = len(environmental_variables)*len(timeperiod)
clim_lists = []
site_lists = []
which_climate_variable = []

for climvar in range(0, len(climate_variables)):
    climate_variable = []
    for tp in range(0, len(timeperiod)):
        for num in range(0, len(sitesARR)):
            climate_file = os.path.join(climate_dir, ('site_'+sitesARR[num]+'_1x1'))
            climate_file = os.path.join(climate_file, (timeperiod[tp][0] + '-climate.nc'))
            clim_lists.append(file_to_mean_standdev(climate_file, climate_variables[climvar], 10))
    which_climate_variable.append(clim_lists)
    clim_lists = []
pr_sq = []
nep_list = []
rh_list = []
npp_list = []
for climvar in range(0, len(climate_variables)):
    pdf = PdfPages('C:/Users/Public/Anaconda/dalton_transect_sites/NEP_'+climate_variables[climvar]+'_ERRORBARS.pdf')
    for tp in range(0, len(timeperiod)):
        for comm in range(0, len(community_type)):
            for env in range(0, len(environmental_variables)):
                sites = []
                if env == 1:
                    nep_site_list = []
                    rh_site_list = []
                    npp_site_list = []
                for num in range(0, len(sitesARR)):
                    site_file = os.path.join(env_dir, ('DHS'+sitesARR[num]+'_'+community_type[comm]))
                    site_file = os.path.join(site_file, (environmental_variables[env]+'_monthly_'+timeperiod[tp][1]+'.nc'))
                    sites.append(file_to_mean_yr(site_file, environmental_variables[env], 10))
                    if env == 1:
                        #print site_lists[comm*climvar*tp + comm*tp + comm]
                        nep_mn = np.subtract(site_lists[comm*climvar*tp + tp*comm+comm][num][0], sites[num][0])
                        nep_yr = np.subtract(site_lists[comm*climvar*tp + tp*comm+comm][num][1], sites[num][1])
                        nep_sd = standard_deviation_multipiece(nep_yr, nep_mn)
                        nep_site_list.append(nep_mn)
                        #print nep_mn
                        rh_site_list.append(sites[num][0])
                       # print sites[num][0]
                        npp_site_list.append(site_lists[comm*climvar*tp + tp*comm+comm][num][0])
                      #  print site_lists[comm*climvar*tp + tp*comm+comm][num][0]
                        fig, ax = plt.subplots()
                        for i in range(0, len(nep_mn)):
                            ax.errorbar(which_climate_variable[climvar][tp*num+num][0][i], nep_mn[i], nep_sd[i], which_climate_variable[climvar][tp*num+num][1][i], fmt='-o', label=str(2010+10*i))
                        #print climvar
                        #print which_climate_variable[climvar][tp*num+num][0]
                        
                        x = np.asarray(which_climate_variable[climvar][tp*num+num][0])
                        y = np.asarray(nep_mn)
                        #regr = LinearRegression(fit_intercept=True)
                        
                        m, b, r, p, s = stats.linregress(x, y)
                        pr_sq.append((r*r, p, ("DHS"+sitesARR[num]+community_type[comm]), climate_variables[climvar]))
                        
                        rise = np.multiply(x, m)
                        rise = np.add(rise, b)
                        ax.plot(x, rise, label='y = '+ str(m)+'x' + str(b))
                        ax.set_xlabel('Annual'+ climate_variables[climvar])
                        ax.set_ylabel('Net Ecosystem Production '+ sites[num][2])
                        ax.set_title('Dalton Highway Site ' + sitesARR[num] +' '+ timeperiod[tp][0]+ ' '+ community_type[comm] + ' 2010-2100')
                        ax.legend(bbox_to_anchor=(.75, .75))
                        pdf.savefig()
                site_lists.append(sites)
            nep_list.append(nep_site_list)
            rh_list.append(rh_site_list)
            npp_list.append(npp_site_list)
    pdf.close()
    '''
    IDEA!!!!!!!!!! SUM ACROSS COMMUNITY TYPES
    '''
diff1 = 0
diff2 = 0
diff3 = 0
diffnpp1 = 0
diffnpp2 = 0
diffnpp3 = 0
diffrh1 = 0
diffrh2 = 0
diffrh3 = 0
sum1 = 0
sum2 = 0
sum3 = 0
for num in range(0, len(sitesARR)):
    site_diff = 0
    site_diff_npp = 0
    site_diff_rh = 0
    site_sum = 0
    site_sum_npp = 0
    site_sum_rh = 0
    for comm in range(0, len(community_type)):
            diff = nep_list[comm][num][len(nep_list[comm][num])-1] - nep_list[comm][num][0]
            #print diff
            diffnpp = npp_list[comm][num][len(npp_list[comm][num])-1] - npp_list[comm][num][0]
            diffrh = rh_list[comm][num][len(rh_list[comm][num])-1] - rh_list[comm][num][0]
            #print int(diffnpp - diffrh) == int(diff)
            #print diff
            site_diff += diff
            site_diff_npp +=diffnpp
            site_diff_rh +=diffrh
            #print site_diff
            fract = diff/nep_list[comm][num][0]
            site_sum += nep_list[comm][num][0]
            percent = fract*100
            if comm == 0:
                diff1 +=diff
                diffnpp1 +=diffnpp
                diffrh1 +=diffrh
                sum1 = nep_list[comm][num][0]
            elif comm == 1:
                diff2 +=diff
                diffnpp2 +=diffnpp
                diffrh2 +=diffrh
                sum2 = nep_list[comm][num][0]
            else:
                diff3 +=diff
                diffnpp3 +=diffnpp
                diffrh3 +=diffrh
                sum3 = nep_list[comm][num][0]
            #print percent
            #print 'Community Type: '+community_type[comm]+' Site: '+sitesARR[num]+' Percent Change: '+str(percent)+'%'
    site_per = (100*site_diff)/site_sum
   #print site_diff
    print ' Site: '+sitesARR[num]+ ' NEP Difference: '+str(site_diff)+ ' Percent: '+str(site_per)+'%'
    print 'NPP Difference: '+str(site_diff_npp)
    print 'RH Difference: ' + str(site_diff_rh)
print 'Difference in NEP CMT04: '+str(diff1)
print 'Difference in NPP CMT04: '+str(diffnpp1)
print 'Difference in RH CMT04: '+str(diffrh1)
print 'Difference in NEP CMT05: '+str(diff2)
print 'Difference in NPP CMT05: '+str(diffnpp2)
print 'Difference in RH CMT05: '+str(diffrh2)
print 'Difference in NEP CMT06: '+str(diff3)
print 'Difference in NPP CMT06: '+str(diffnpp3)
print 'Difference in RH CMT06: '+str(diffrh3)
        
dtype = [('r_sq', float),('p', float),('where', 'S10'),('climate', 'S10')]
values = pr_sq
pr_sq = np.array(values, dtype=dtype)

print np.sort(pr_sq, order='r_sq')
