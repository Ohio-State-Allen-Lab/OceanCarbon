#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 13:27:19 2021

@author: AbbieEnders
"""
# adjusting the 
import pandas as pd
import os
import math
import numpy
import re
import glob
# rotating matrix
def rotate_180(array, M, N, out):
    for i in range(M):
        for j in range(N):
            out[i, N-1-j] = array[M-1-i, j]

# equation
# Ci = GCz(Cp/(Ck,inges + Cp))(1-a)(pi%)(ti)(o)
# variables
# Ci = carbon of ith species
# G  = 1/day zooplanktonic growth rate
# Cz = uM carbon zooplankton concentration
# Cp = uM carbon planktonic carbon atom concentration
# Ck,inges = uM carbon half saturation for ingestion
# a = assimilation efficiency
# pi% = percentage of the ith macromolecule content in a typical planktonic cell
# ti = day lifetime of the ith macromolecule, total restricted to 2
# o = coating of the surface based on partial adsorption
# ChlA = remotely sensed by NASA MODIS
# Cmr = chlorophyll mass ratio (multiply ChlA by 50 to get Cp)

# variables defined
a = 0.75
G = 1
#Cz = 0.5
Ckinges = 7
tprot = 10
tlip = 2
pprot = 0.6
plip = 0.2
Cmr = 50 
C_ratio = 0.5
CpRef = 10 # carbon protein reference uM
ClRef = 0.5 # carbon lipid reference uM
np = 0.5
nl = 1
ap = 1
al = 1


mol_mass_carbon = 12.01 # g/mol
ocean_surf_area = 3.60580510*10**14  #m^2

earth_surf_area = 5.10082000*10**14 # m^2
num_of_instances = 180 *360 #1 steps in lat and lon
a_pixel = earth_surf_area/num_of_instances
surfprot = 0.002*a_pixel # grams
surflip = 0.0025*a_pixel # grams

# relate num of instances over ocean surface area

# new dataframe
# math
# if ChlA = 99999 do nothing
path = '/Users/AbbieEnders/Library/CloudStorage/OneDrive-TheOhioStateUniversity/nasa modis data'
os.chdir(path)
all_files_chl = glob.glob(path + "/*chl.csv")
#print(all_files_chl)
all_files_zoo = glob.glob(path + "/*zoo.csv")
#lat_list = list(range(-90, 91,0.5))
lat_list = numpy.arange(-90,91,0.5).tolist()
#lat_list.append(list(range(90,-1)))
#print(lat_list)
os.chdir(path)
list_sum = []
for filename in all_files_chl:
	
    data = pd.read_csv(filename, index_col=0)
    Cp_temp = Cmr*data
    Cz = re.sub('chl','zoo',filename)
    zoo = pd.read_csv(Cz, index_col = 0)
    zoo = zoo/1000
    C_prot_temp = G*zoo*(Cp_temp/(Ckinges + Cp_temp))*(1-a)*(pprot)*(tprot)
    C_lip_temp = G*zoo*(Cp_temp/(Ckinges + Cp_temp))*(1-a)*(plip)*(tlip)
    # nfn = re.sub('chl','prot',filename)
    # C_prot_temp.to_csv(nfn)
    C_sum_temp = C_prot_temp + C_lip_temp
    c_lip_sum = C_sum_temp.iloc[20:25,160:170]
    c_lip_sum = c_lip_sum.sum()
    c_lip_sum = c_lip_sum.sum()
    list_sum.append(c_lip_sum)
    #print(filename, c_lip_sum)
    theta_prot = (((1/CpRef)**np)*((ap*C_prot_temp)**np))/(1 +((((1/CpRef)**np)*((ap*C_prot_temp)**np))+(((1/ClRef)**nl)*((al*C_lip_temp)**nl))))
    theta_lip =  (((1/ClRef)**nl)*((al*C_lip_temp)**nl))/(1 +((((1/CpRef)**np)*((ap*C_prot_temp)**np))+(((1/ClRef)**nl)*((al*C_lip_temp)**nl))))
    theta = theta_prot + theta_lip
    # nfn = re.sub('chl','theta',filename)
    # theta.to_csv(nfn)
    sums = 0
    counter = 0
    my_dict = {}
    for i in lat_list:
        i=float(i)
        C_total = theta_lip.iloc[counter,:]*surflip*math.cos(numpy.deg2rad(abs(i))) + theta_prot.iloc[counter,:]*surfprot*math.cos(numpy.deg2rad(abs(i)))
        my_dict[i]=(abs(C_total))
        sums = sums + C_total.to_numpy().sum()
        counter += 1
	# get sum of carbon in region with measurements
   
    # savefile = re.sub('chl','carbonn',filename)	
    # C_total = pd.DataFrame.from_dict(my_dict, orient = 'index')
    # C_total.to_csv(savefile)
	#print(sums)
	#C_total = C_total.fillna(0)
	#sum_total_carbon = C_total.to_numpy().sum()
	#print(sum_total_carbon)
	#Gt_c = sum_total_carbon * (1*10**-15) 
	#print(Gt_c, savefile)
    
    

total = sum(list_sum)/12
print(total)
#data_feb = pd.read_csv('/Users/AbbieEnders/Library/CloudStorage/OneDrive-TheOhioStateUniversity/nasa modis data/feb_2005_chl.csv',index_col=0)


