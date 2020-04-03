#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr  1 16:09:32 2020

@author: samantha
"""

# # Homework 7 Template
# 
# Rixin Li & G . Besla
# 
# look for "****"



# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# matplotlib provides powerful functions for plotting figures
import matplotlib.pyplot as plt
# astropy provides unit system and constants for astronomical calculations
import astropy.units as u
import astropy.constants as const
# import Latex module so we can display the results with symbols
from IPython.display import Latex

# **** import CenterOfMass to determine the COM pos/vel of M33
from CenterOfMass import CenterOfMass

# **** import the GalaxyMass to determine the mass of M31 for each component
from GalaxyMass import ComponentMass

# # M33AnalyticOrbit

# In[ ]:


class M33AnalyticOrbit:
    """ Calculate the analytical orbit of M33 around M31 """
    
    def __init__(self, filename): 
        """inputs:
            filename for file in which we will store the integrated orbit""" 

        ### get the gravitational constant (the value is 4.498502151575286e-06)
        self.G = const.G.to(u.kpc**3/u.Msun/u.Gyr**2).value
        
        
        
        ### **** store the output file name
        self.fileout = filename
        
        ### get the current pos/vel of M33 
        # **** create an instance of the  CenterOfMass class for M33 
        
        M33COM = CenterOfMass("M33_000.txt", 2)
        #print("got M33COM")
        # **** store the position VECTOR of the M33 COM (.value to get rid of units)
        M33_COMP = M33COM.COM_P(0.1)
        
        # **** store the velocity VECTOR of the M33 COM (.value to get rid of units)
        M33_COMV = M33COM.COM_V(M33_COMP[0], M33_COMP[1], M33_COMP[2])
        
        ### get the current pos/vel of M31 
        # **** create an instance of the  CenterOfMass class for M31 
        M31COM = CenterOfMass("M31_000.txt", 2)
        #print("got M31COM")
        # **** store the position VECTOR of the M31 COM (.value to get rid of units)
        M31_COMP = M31COM.COM_P(0.1)
        # **** store the velocity VECTOR of the M31 COM (.value to get rid of units)
        M31_COMV = M31COM.COM_V(M31_COMP[0], M31_COMP[1], M31_COMP[2])
        
        ### store the DIFFERENCE between the vectors posM33 - posM31
        # **** create two VECTORs self.r0 and self.v0 and have them be the
        # relative position and velocity VECTORS of M33
        #M31_COMP is a vector
        self.r0 = (M33_COMP - M31_COMP).value
        self.v0 = (M33_COMV - M31_COMV).value
        
        ### get the mass of each component in M31 
        ### disk
        # **** self.rdisk = scale length (no units)
        self.rdisk = 5 #in kpc
        
        
        # **** self.Mdisk set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mdisk = ComponentMass("M31_000.txt", 2)*1e12 #in M_sun
        ### bulge
        # **** self.rbulge = set scale length (no units)
        self.rbulge = 1
        # **** self.Mbulge  set with ComponentMass function. Remember to *1e12 to get the right units Use the right ptype
        self.Mbulge = ComponentMass("M31_000.txt", 3)*1e12
        # Halo
        # **** self.rhalo = set scale length from HW5 (no units)
        self.rhalo = 62 #kpc
        # **** self.Mhalo set with ComponentMass function. Remember to *1e12 to get the right units. Use the right ptype
        self.Mhalo = ComponentMass("M31_000.txt", 1)*1e12
        
    
    def HernquistAccel(self, M, r_a, vector): # it is easiest if you take as an input the position VECTOR 
        """ Inputs:
                M = total halo or bulge mass
                r_a = corresponding scale length
                vector = position vector
            Output: acceleration vector for Hernquist potential"""
        
        ### **** Store the magnitude of the position vector
        rmag = comp_mag(vector)
        #print(vector)
        ### *** Store the Acceleration
        Hern =  -((self.G*M) / (rmag*(r_a + rmag)**2)) * vector
        # NOTE: we want an acceleration VECTOR so you need to make sure that in the Hernquist equation you 
        # use  -G*M/(rmag *(ra + rmag)**2) * r --> where the last r is a VECTOR 
        
        return Hern
    
    
    
    def MiyamotoNagaiAccel(self, M, rd, vector):# it is easiest if you take as an input a position VECTOR  r 
        """ approximation that mimics the exponential disk profile at 
        distances far from the disk.
        Inputs: 
            M = mass 
            rd = disk radius of M31
            vector = position vector
        Output:
            Acceleration vector"""
    
        
        ### Acceleration **** follow the formula in the HW instructions
        # AGAIN note that we want a VECTOR to be returned  (see Hernquist instructions)
        # this can be tricky given that the z component is different than in the x or y directions. 
        # we can deal with this by multiplying the whole thing by an extra array that accounts for the 
        # differences in the z direction:
        #  multiply the whle thing by :   np.array([1,1,ZSTUFF]) 
        # where ZSTUFF are the terms associated with the z direction
        zd = self.rdisk/5.0
        R = np.sqrt(vector[0]**2 + vector[1]**2)
        B = rd + np.sqrt(vector[2]**2 + zd**2)        
        #MNAccel = (-(self.G*M) / (R**2 + B**2)**(1.5)) * np.array([1, 1, B/(vector[2]**2 + zd**2)**(1/2)])
        return ((-self.G * M)/(R**2 + B**2)**1.5)*vector*np.array([1, 1, B/np.sqrt(vector[2]**2 + zd**2)])
        # the np.array allows for a different value for the z component of the acceleration
     
    
    def M31Accel(self, vector): # input should include the position vector, r
        """ Inputs:
            vector = 3D position vector
        Output: 
            3D vector of the total acceleration"""

        ### Call the previous functions for the halo, bulge and disk
        # **** these functions will take as inputs variable we defined in the initialization of the class like 
        # self.rdisk etc. 
        halo_accel = self.HernquistAccel(self.Mhalo, self.rhalo, vector)
        bulge_accel = self.HernquistAccel(self.Mbulge, self.rbulge, vector)
        disk_accel = self.MiyamotoNagaiAccel(self.Mdisk, self.rdisk, vector)
            # return the SUM of the output of the acceleration functions - this will return a VECTOR 
        return np.sum([halo_accel, bulge_accel, disk_accel], axis = 0)
    
    
    
    def LeapFrog(self, dt, r, v): # take as input r and v, which are VECTORS. Assume it is ONE vector at a time
        """Inputs:
            dt = time interval for integration
            r = position vector
            v = position vecor
            Output:
                new position and velocity vectors"""
        
        # predict the position at the next half timestep
        rhalf = r + v*dt/2
        
        # get acceleration at half step
        a = self.M31Accel(rhalf)
        
        # predict the final velocity at the next timestep using the acceleration field at the rhalf position 
        vnew = v + a*dt
        
        # predict the final position using the average of the current velocity and the final velocity
        # this accounts for the fact that we don't know how the speed changes from the current timestep to the 
        # next, so we approximate it using the average expected speed over the time interval dt. 
        rnew = r + dt*(v + vnew)/2
        
        return rnew, vnew# **** return the new position and velcoity vectors
    
    
    
    def OrbitIntegration(self, t0, dt, tmax):
        """ Input: 
                t0 = initial starting time
                dt = time interval
                tmax = final time
            Output:
                there is no output"""
        #starting COM position and velocity of M33 relative to M31
        r = self.r0
        v = self.v0
        
        # initialize the time to the input starting time
        t = t0 
        
        # initialize an empty array of size :  rows int(tmax/dt)+1  , columns 7
        size = np.int(tmax/dt) + 1
        orbit = np.zeros((size, 7))
        
        # initialize the first row of the orbit
        orbit[0] = t0, *tuple(self.r0), *tuple(self.v0)
        # this above is equivalent to 
        # orbit[0] = t0, self.r0[0], self.r0[1], self.r0[2], self.v0[0], self.v0[1], self.v0[2]
        
        
        # initialize a counter for the orbit.  
        i = 1 # since we already set the 0th values, we start the counter at 1
        
        # start the integration (advancing in time steps and computing LeapFrog at each step)
        while (t < tmax):  # as long as t has not exceeded the maximal time 
            
            # **** advance the time by one timestep, dt
            t += dt
            
            # ***** advance the position and velocity using the LeapFrog scheme
            # remember that LeapFrog returns a position vector and a velocity vector  
            # as an example, if a function returns three vectors you would call the function and store 
            # the variable like:     a,b,c = function(input)
            r, v = self.LeapFrog(dt, r, v)
         
            # **** store the new information in the first column of the ith row
            orbit[i] = t, *tuple(r), *tuple(v)

            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            # TIP:  if you want columns 5-7 of the Nth row of an array called A, you would write : 
            # A[n, 5:8] 
            # where the syntax is row n, start at column 5 and end BEFORE column 8
            
            
            # ****  store the new position vector into the columns with indexes 1,2,3 of the ith row of orbit
            
            
            # **** update counter i , where i is keeping track of the number of rows (i.e. the number of time steps)
            i += 1
        
        
        # write the data to a file
        np.savetxt(self.fileout, orbit, fmt = "%11.3f"*7, comments='#', 
                   header="{:>10s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}{:>11s}"\
                   .format('t', 'x', 'y', 'z', 'vx', 'vy', 'vz'))
        
        # there is no return function
        
# Plot
    def ReadOrbits(self):
        """ Read in file that was output so we can plot """
        return np.genfromtxt(self.fileout, dtype = None, names = True)

def RelativeMag(orbit1, orbit2):
    """ Function to compute mag of difference btw two vectors
    Inputs:
        orbit1, orbit2 = 3D position and velocity vectors
    Output:
        magnitude of relative position and velocities"""
        
    r = np.sqrt((orbit1['x'] - orbit2['x'])**2 + (orbit1['y'] - orbit2['y'])**2 + (orbit1['z'] - orbit2['z'])**2)
    v = np.sqrt((orbit1['vx'] - orbit2['vx'])**2 + (orbit1['vy'] - orbit2['vy'])**2 + (orbit1['vz'] - orbit2['vz'])**2)
    return r, v
            
def comp_mag(vector):
    return np.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2)     


OrbitOut = "orbitout.txt"

# create an instance of the orbit
M33_Accel = M33AnalyticOrbit(OrbitOut)

M33_Accel_Vector = M33_Accel.M31Accel(M33_Accel.r0)

M33_Accel.OrbitIntegration(0, 0.05, 10)
# t0 = 0 Gyr
# dt = 0.05 Gyr
# tmax = 10 Gyr

# file needs to be read in and put into an array
M33_Orbits = M33_Accel.ReadOrbits()
mag_position_orbits = comp_mag([M33_Orbits['x'], M33_Orbits['y'], M33_Orbits['y']])
mag_velocity_orbits = comp_mag([M33_Orbits['vx'], M33_Orbits['vy'], M33_Orbits['vz']])      

# Overplot HW6 assignment for M33's orbit wrt M31 from the simulation
M31_Orbit = np.genfromtxt("Orbit_M31.txt", dtype = None, names = True)
M33_Orbit = np.genfromtxt("Orbit_M33.txt", dtype = None, names = True)
            
            
# determine magnitude of relative position and velocities
M33_M31_position, M33_M31_velocity = RelativeMag(M31_Orbit, M33_Orbit)


# Graphs: Position
# Predicted M33 orbit
fig = plt.figure(figsize = (5, 5))
ax = plt.subplot(111)

plt.plot(M33_Orbits['t'], mag_position_orbits, color = 'blue', linewidth = 5, label = 'M33-M31 Analytical') 

plt.xlabel('Time (Gyr)')
plt.ylabel('Separation (kpc)')

legend = ax.legend(loc = 'upper left', fontsize = 'small')

plt.title('M33 motion relative to M31', fontsize = 22)

# overplot M33's orbit wrt M31
ax.plot(M31_Orbit['t'], M33_M31_position, color = 'red', linewidth = 5, label = 'M333-M31 Simulation')
legend = ax.legend(loc = 'upper left', fontsize = 'small')
plt.show
plt.savefig('M33 position relative to M31', dpi = 350)
            
# Graphs: Velocity
# Predicted M33 orbit
fig = plt.figure(figsize = (5, 5))
ax = plt.subplot(111)

plt.plot(M33_Orbits['t'], mag_velocity_orbits, color = 'blue', linewidth = 5, label = 'M33-M31 Analytical') 

plt.xlabel('Time (Gyr)')
plt.ylabel('Velocity (km/s)')

legend = ax.legend(loc = 'upper left', fontsize = 'small')

plt.title('M33 velocity relative to M31', fontsize = 22)

# overplot M33's orbit wrt M31
ax.plot(M31_Orbit['t'], M33_M31_velocity, color = 'red', linewidth = 5, label = 'M333-M31 Simulation')
legend = ax.legend(loc = 'upper left', fontsize = 'small')
plt.show
plt.savefig('M33 velocity relative to M31', dpi = 350)
            
""" QUESTIONS:
1. How do the plots compare?
    Both the velocity and position plots are very similar for the first
    few Gyr. However, the simulation begins to oscillate much more than
    the analytical. The analytical is simply one large oscillation. 
    
2. What missing physics could make the difference?
    The physics mising could be that the analytical solution does not
    include the presence of the MW. Another large galaxy, similar to M31
    would greatly impact the position and velocity of M31 and M33. 
    
3. How might you include the effects of the MW?
    The MW could be included as another body, and the calculations redone
    with the COM being that of the system instead of instead of just M31 and
    M33. We would also need to account for the acceleration due to MW's 
    gravity.
            
            
            
            
            
            
            