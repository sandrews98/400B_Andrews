# Homework 4
# Center of Mass Position and Velocity
# Samantha Andrews
#2/9/20

###############################################################################
# Keep in mind this is just a template; you don't need to follow every step and feel free to change anything.
# We also strongly encourage you to try to develop your own method to solve the homework.
###############################################################################

# import modules
import numpy as np
import astropy.units as u
from ReadFile import Read
import math


class CenterOfMass:
# Class to define COM position and velocity properties of a given galaxy 
# and simulation snapshot
    
    
    def __init__(self, filename, ptype):
    # Initialize the instance of this Class with the following properties:
    
        # read data in the given file using Read
        self.time, self.total, self.data = Read(filename)                                                                                             

        #create an array to store indexes of particles of desired Ptype                                
        self.index = np.where(self.data['type'] == ptype)
        #print("Index = ", self.index)

        # store the mass, positions, velocities of only the particles of the given type
        # the following only gives the example of storing the mass
        self.m = self.data['m'][self.index]
        self.x = self.data['x'][self.index]
        #print("x[0-4]",self.x[0],self.x[1],self.x[2],self.x[3],self.x[4])
#        print(self.x)
#        print("Length of x",len(self.x))
        self.y = self.data['y'][self.index]
        self.z = self.data['z'][self.index]
        self.vx = self.data['vx'][self.index]
        self.vy = self.data['vy'][self.index]
        self.vz = self.data['vz'][self.index]
        #self refers to quantities that are common to the object
        #each function created must start with self as an input
        
    def comp_mag(self,x,y,z):
        #fucntion to compute magnitude given the coordinates (x,y,z) in kpc
        #   from the center of the mass position
        # NONE OF THESE VALUES ARE GLOBAL FOR THE CLASS SO I AM NOT USING SELF IN FRONT OF THEM   
        #initialize
        mag=0
        #Compute mag
        mag=math.sqrt(x**2+y**2+z**2)
        #Returen the mag
        return mag
        
        
    def COMdefine(self,a,b,c,m):
    # Function to compute the center of mass position or velocity generically
    # input: array (a,b,c) of positions or velocities and the mass
    # returns: 3 floats  (the center of mass coordinates)
    #values not global for the class so not using self in front of them

        # write your own code to compute the generic COM using Eq. 1 in the homework instructions
        # xcomponent Center of mass
        Acom = np.sum(np.multiply(a,m))/np.sum(m)
        # ycomponent Center of mass
        Bcom = np.sum(np.multiply(b,m))/np.sum(m)
        # zcomponent Center of mass
        Ccom = np.sum(np.multiply(c,m))/np.sum(m)
        #print("Acom", Acom, "Bcom", Bcom, "Ccom", Ccom)
        return Acom, Bcom, Ccom
    
    def COM_P(self, delta, VolDec):
    # Function to specifically return the center of mass position and velocity                                         
    # input:                                                                                                           
    #        particle type (1,2,3)                                                                                     
    #        delta (tolerance) decides whether the COM position has converged                                                                                        
    # returns: One vector, with rows indicating:                                                                                                                                                                            
    #       3D coordinates of the center of mass position (kpc)                                                             

        # Center of Mass Position                                                                                      
        ###########################                                                                                    

        # Try a first guess at the COM position by calling COMdefine
        #goal is to refine the COM position calc iteratively to make sure the position has converged
        XCOM, YCOM, ZCOM = self.COMdefine(self.x, self.y, self.z, self.m)
        # compute the magnitude of the COM position vector.
        #when referring to functions already defined in class use self.function
        RCOM = self.comp_mag(XCOM, YCOM, ZCOM)
        #print("First RCOM", RCOM)


        # iterative process to determine the center of mass                                                            

        # change reference frame to COM frame                                                                          
        # compute the difference between particle coordinates                                                          
        # and the first guess at COM position
        # subtract first guess COM position from partivle position
        xNew = self.x - XCOM
        yNew = self.y - YCOM
        zNew = self.z - ZCOM
        
        #create an array to store the mag of the new position vecors of all
        #   particles in the COM frame
        RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)
        
        #sanity check
#        print("xNew", xNew)
#        print("yNew", yNew)
#        print("zNew", zNew)
#        print("first value", np.sqrt(xNew[0]**2 + yNew[0]**2 + zNew[0]**2))
#        print("RNEW", RNEW)

        # find the max 3D distance of all particles from the guessed COM                                               
        # will re-start at half that radius (reduced radius)                                                           
        RMAX = max(RNEW)/VolDec
        #print(f"RMax = {RMAX}")
        
        # pick an initial value for the change in COM position                                                      
        # between the first guess above and the new one computed from half that volume
        # it should be larger than the input tolerance (delta) initially
        CHANGE = 1000.0 #[kpc]

        # start iterative process to determine center of mass position                                                 
        # delta is the tolerance for the difference in the old COM and the new one.    
        
        while (CHANGE > delta): #delta is the tolerance
            # select all particles within the reduced radius (starting from original x,y,z, m)
            #runs as long as the difference btw RCOM and a new COM position
            #   is larger than some tolerance (delta)
            # write your own code below (hints, use np.where)
            index2 = np.where(RNEW < RMAX)
            """
            #show that indexing is narrowing in sync
            print("New loop through while")
            print("RNew = ", RNEW)
            print("RMax = ", RMAX)
            print("index2 = ", index2)
            print("RNEW[0-4] = ", RNEW[index2][0], RNEW[index2][1], RNEW[index2][2], RNEW[index2][3], RNEW[index2][4])
            """
            x2 = self.x[index2]
            y2 = self.y[index2]
            z2 = self.z[index2]
            m2 = self.m[index2]
            
            #Check values are narrowing
            """
            print("Len x2", len(x2))
            print("x2[0-4]", x2[0], x2[1], x2[2], x2[3], x2[4])
            print("y2[0-4]", y2[0], y2[1], y2[2], y2[3], y2[4])
            print("z2[0-4]", z2[0], z2[1], z2[2], z2[3], z2[4])
            print("m2[0-4]", m2[0], m2[1], m2[2], m2[3], m2[4])
            """
            # Refined COM position:                                                                                    
            # compute the center of mass position using                                                                
            # the particles in the reduced radius
            # write your own code below
            XCOM2, YCOM2, ZCOM2 = self.COMdefine(x2, y2, z2, m2)
            #print("Xcom2", XCOM2, "Ycom2", YCOM2, "Zcom2", ZCOM2)
            
            # compute the new 3D COM position
            # write your own code below
            RCOM2 = self.comp_mag(XCOM2, YCOM2, ZCOM2)
            #print(f"new COM magnitude = {RCOM2}")

            # determine the difference between the previous center of mass position                                    
            # and the new one.                                                                                         
            CHANGE = np.abs(RCOM - RCOM2)
            # uncomment the following line if you wnat to check this                                                                                               
            #print ("CHANGE = ", CHANGE)                                                                                     

            # Before loop continues, reset : RMAX, particle separations and COM                                        

            # reduce the volume by a factor of 2 again                                                                 
            RMAX = RMAX/VolDec
            # check this.                                                                                              
            #print ("maxR", maxR)                                                                                      

            # Change the frame of reference to the newly computed COM.                                                 
            # subtract the new COM
            # write your own code below
            xNew = self.x - XCOM2
            yNew = self.y - YCOM2
            zNew = self.z - ZCOM2
            RNEW = np.sqrt(xNew**2 + yNew**2 + zNew**2)
            #new array of magnitudes for new from of reference COM
            #print("RNEW", RNEW)

            #refine volume again
            #RMAX = RMAX/2.
            #print("RMAX", RMAX)

            # set the center of mass positions to the refined values                                                   
            XCOM = XCOM2
            YCOM = YCOM2
            ZCOM = ZCOM2
            RCOM = RCOM2

            # create a vector to store the COM position                                                                                                                                                       
            #COMP = [XCOM, YCOM, ZCOM]

        # set the correct units usint astropy and round all values
        # and then return the COM positon vector
        # write your own code below
        COMP = np.around([XCOM, YCOM, ZCOM], 2)*u.kpc
        return COMP
    

    def COM_V(self, COMX,COMY,COMZ):
        # Center of Mass velocity
        # input: X, Y, Z positions of the COM
        # returns 3D Vector of COM Velocities
        
        # the max distance from the center that we will use to determine the center of mass velocity                   
        RVMAX = 15.0*u.kpc

        # determine the position of all particles relative to the center of mass position
        # write your own code below
        xV = self.x - COMX
        yV = self.y - COMY
        zV = self.z - COMZ
        RV = np.sqrt(xV**2 + yV**2 + zV**2)
        
        # determine the index for those particles within the max radius
        # write your own code below
        indexV = np.where(RV < RVMAX.value)

        # determine the velocity and mass of those particles within the mas radius
        # write your own code below
        vxnew = self.vx[indexV]
        vynew = self.vy[indexV]
        vznew = self.vz[indexV]
        mnew = self.m[indexV]
        
        # compute the center of mass velocity using those particles
        # write your own code below
        VXCOM, VYCOM, VZCOM = self.COMdefine(vxnew, vynew, vznew, mnew)

        # create a vector to store the COM velocity
        # set the correct units usint astropy
        # round all values
        # write your own code below
        COMV = np.around([VXCOM, VYCOM, VZCOM], 2)*u.km/u.s

        # return the COM vector                                                                                        
        return COMV
    
"""

# ANSWERING QUESTIONS
#######################


# Create  a Center of mass object for the MW, M31 and M33
# below is an example of using the class for MW
MWCOM = CenterOfMass("MW_000.txt", 2)
# now write your own code to answer questions
#QUESTION 1:
# below gives you an example of calling the class's functions
# MW:   store the position and velocity COM 
MW_COMP = MWCOM.COM_P(0.1)
print("MW_COMP", MW_COMP)

MW_COMV = MWCOM.COM_V(MW_COMP[0].value,MW_COMP[1].value,MW_COMP[2].value)
print("MW_COMV", MW_COMV)


#for M31
M31COM = CenterOfMass("M31_000.txt", 2)
M31_COMP = M31COM.COM_P(0.1)
print("M31_COMP", M31_COMP)

M31_COMV = M31COM.COM_V(M31_COMP[0].value, M31_COMP[1].value, M31_COMP[2].value)
print("M31_COMV", M31_COMV)

#for M33
M33COM = CenterOfMass("M33_000.txt", 2)
M33_COMP = M33COM.COM_P(0.1)
print("M33_COMP", M33_COMP)

M33_COMV = M33COM.COM_V(M33_COMP[0].value, M33_COMP[1].value, M33_COMP[2].value)
print("M33_COMV", M33_COMV)

#QUESTION 2:
#What is the mag of the current separation and velocity between MW and M31?
sep = MW_COMP - M31_COMP
mag_sep = np.sqrt(sep[0]**2 + sep[1]**2 + sep[2]**2)
print(f"The magnitude fo the separation between MW and M31 is {mag_sep}")

velocity = MW_COMV - M31_COMV
mag_velocity = np.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2)
print(f"The magnitude of the velocity between MW and M31 9s {mag_velocity}")

#QUESTION 3:
#What is the magnitude of the current separation and velocity between M33 and M31?
sep = M33_COMP - M31_COMP
mag_sep = np.sqrt(sep[0]**2 + sep[1]**2 + sep[2]**2)
print(f"The magnitude fo the separation between MW and M31 is {mag_sep}")

velocity = M33_COMV - M31_COMV
mag_velocity = np.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2)
print(f"The magnitude of the velocity between MW and M31 9s {mag_velocity}")

#QUESTION 4:
#It is important to use the iterative process for this system because
#   MW and M31 are about to collide so there will be a lot of changes
#   to the distribution of the mass. It is important to have the disk movements 
#   since the halo will move around. We need to know the center of the stars, 
#   and the iterative process allows us to do that. 
"""
