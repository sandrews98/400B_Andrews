#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 12:25:38 2020

@author: samantha
"""

"""
This is the code for HW2 in ASTR400B
It will read in the MW_000.txt data file
"""

import numpy as np
import astropy.units as u

#define a function to read in data file
def Read(filename):
    file = open(filename, 'r') #opens file to read
    
    #read in first line and store w/ units of 10Myr
    line1 = file.readline()
    label, value = line1.split()
    time = float(value)*u.Myr
    
    #now read in second line for total number of particles
    line2 = file.readline()
    label2, number_particles = line2.split()
    
    #close file
    file.close()
    
    #store rest of file
    data = np.genfromtxt(filename, dtype = None, names = True, skip_header = 3)
    
    """
    dtype = None means line is split using the white spaces
    skip_header = skips the first three lines
    name = True creates arrays to store the data with labels "type, m, x, y, z, 
    vx, vy, vz
    """

    return time, number_particles, data

#test
def TestRead():
    time = 0
    particles = 0
    data = []
    filename = 'MW_000.txt'
    time, particles, data = Read(filename)
    print(f"Time = {time}\n Number of particles = {particles}")
    print(f"first type = {data['type'][0]} mass = {data['m'][0]} x = {data['x'][0]}")
    
    return
    
TestRead()
    