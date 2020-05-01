
# coding: utf-8

# In[ ]:


# In Class Lab 6
# Surface Brightness Profiles


# In[ ]:


# Load Modules
import numpy as np
import astropy.units as u

# import plotting modules
import matplotlib.pyplot as plt
import matplotlib

# my modules
from ReadFile import Read
from CenterOfMass import CenterOfMass
from MassProfile import MassProfile


# # Lab 6: Sersic Profiles
# 
# In this lab we will use Homework 5 solutions to compute the mass profile of M31's bulge. 
# We will turn the mass profile into a density profile and see if we can fit it reasonably well 
# with a sersic profile. 

# # Part A : 
# 
# Create a function called `SersicE` that returns the Sersic Profile in terms of the effective radius $R_e$ (i.e. the half light radius).
# 
# $I(r) = I_e exp^{-7.67 ( (r/R_e)^{1/n} - 1)}$
# 
# Where 
# 
# $ L = 7.2 I_e \pi R_e^2$
# 
# and  $R_e$ is the half light radius.  We will assume a mass to light ratio for the stellar bulge of 1, so this is also the half mass radius.
# 
# The function should take as input: the radius, $R_e$, $n$ and the total stellar mass of the system.
# 

# In[ ]:
#give an array of radii
#function to return Sersic profile for an elliptical system
def SersicE(R, Re, n, Mtot):
    #Inputs:
    #   R: array of radii [kpc]
    #   Re: half mass radius [kpc]
    #   n: Sersic index
    #   Mtot, total stellar mass [Msun]
    #Return:
    #   surface brightness profile [Lsun/kpc^2]
    #asuume M/L = 1
    L = Mtot #total luminosity = total stellar mass
    #Ie = L/7.2 pi Re^2
    Ie = L/7.2/np.pi/Re**2
    
    #exponent
    A = (R/Re)**(1/n) - 1
    
    return Ie*np.exp(-7.67*A)


  


# # Part B
# 
# Compute the surface density profile for the simulated bulge
# 
# a) Create an instance of the MassProfile Class for M31. Store it as a variable `M31`. 
# 
M31 = MassProfile("M31", 0)

# b) Create an array of radii from 0.1 kpc to 30 kpc in increments of 0.1
# 
R = np.arange(0.1, 30, 0.1)
BulgeMass = M31.MassEnclosed(3, R)  #particle type 3 = bulge
print(BulgeMass[10])

# c) Define a new array called `BulgeMass`, that uses the function `MassEnclosed` within MassProfile to compute the mass profile of the bulge.  Get rid of astropy units in `BulgeMass` by adding `.value` 
# I = L / 4 pi D^2
BulgeI = BulgeMass/4/np.pi/R**2

# d) Compute the surface mass density profile for the simulated bulge and store it as an array called `BulgeI`. Assuming M/L ~ 1 this is also the surface brightness profile in Lsun/kpc^2

# # Part C
# 
# Compute $R_e$, the half mass radius, for the bulge
BulgeTotal = np.max(BulgeMass) #total mass of bulge
Low = BulgeTotal/2
High = BulgeTotal/2 + BulgeTotal/2*0.01
index = np.where((BulgeMass > Low) & (BulgeMass < High))
Re = R[index] #arrays are listed in the same order so you can use same index
print(BulgeTotal)
print(BulgeTotal/2)

# # Part D
# 
# a) Plot the surface density profile of the simulated bulge
# 
# b) Plot the Sersic profile, assuming a de Vaucouleurs Profile.
# 
# c) If the profiles don't match, try changing either $R_e$ or $n$

# In[ ]:


# Plot the Bulge density profile vs 
# the Sersic profile
####################################


fig = plt.figure(figsize=(10,10))
ax = plt.subplot(111)


# plot the bulge luminosity density as a proxy for surface brighntess
plt.semilogy(R,BulgeI, color='black',linewidth=3, label='Bulge Density')


# YOU ADD HERE: Sersic fit to the surface brightness Sersic fit
# Sersic 
plt.semilogy(R, SersicE(R, Re, 5.45, BulgeTotal), color = 'red', linestyle = "-.", Linewidth = 3, label = 'Sersic n = 5.54')


#adjust tick label font size
label_size = 22
matplotlib.rcParams['xtick.labelsize'] = label_size 
matplotlib.rcParams['ytick.labelsize'] = label_size


# Add axis labels
plt.xlabel('Radius (kpc)', fontsize=22)
plt.ylabel('Log(I)  $L_\odot/kpc^2$', fontsize=22)



# add a legend with some customizations.
legend = ax.legend(loc='upper right',fontsize='x-large')

plt.show()
