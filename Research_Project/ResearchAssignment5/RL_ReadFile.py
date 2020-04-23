# import necessary modules
# numpy provides powerful multi-dimensional arrays to hold and manipulate data
import numpy as np
# pandas provides optimal functions to read column data from files
import pandas as pd
# astropy provides unit system for astronomical calculations
import astropy.units as u

"""
Define a function that will read in a file
Input:  filename, e.g. "MW_000.txt"
Returns:  time (in Myr), total number of particles 
          and an array with data
USAGE :   time, total, data = Read("filename")
"""

def Read(filename):
    """ Read particle data from a file """
    
    # open the file with read-only mode
    f = open(filename, 'r')
    
    # read in the headers line by line with readline()
    # for each line, split() will split text by white spaces
    # from the first line, we extract the time info of one snapshot
    # from the second line, we extract the total number of particles
    time = float(f.readline().split()[1]) * u.Myr
    total = int(f.readline().split()[1])
    
    # close the file object
    f.close()
    
    # genfromtxt() reads data line by line, which can be optimized more
    # the code below returns the same thing but is ~10 times faster
    #data = np.genfromtxt(filename,dtype=None,names=True,skip_header=3)
    
    # read_csv() from pandas library can read column data super fast
    # here we provide filename, and specify whitespace as the delimeter
    # also, we skip the first four rows and will give columns names later;
    # read_csv() returns a panda DataFrame, but we can use ".values" to obtain
    # the internal numpy array, and then np.ascontiguousarray() will re-organize
    # the numpy array in rows (was in columns) which will be helpful later
    data = np.ascontiguousarray(pd.read_csv(filename, delim_whitespace=True, 
                                            skiprows=4, header=None).values)
    
    # now that data is stored in memory row by row (particle by particle),
    # we are able to assign a new data type to this data object,
    # which essentially groups one row (8 numbers) into one type structure;
    # in this way, we are also able to use the column names in header, 
    # like: data['m']. There are multiple ways to achieve this,
    # the first two are equivalent, the last one groups position and velocity
    
    #data.dtype = np.dtype([('type', 'f8'), ('m', 'f8'), 
    #                       ('x', 'f8'), ('y', 'f8'), ('z', 'f8'), 
    #                       ('vx', 'f8'), ('vy', 'f8'), ('vz', 'f8')])
    #data.dtype = {'names': ('type', 'm', 'x', 'y', 'z', 'vx', 'vy', 'vz'), 
    #              'formats': ('f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8', 'f8')}
    data.dtype = np.dtype([('type', 'f8'), ('m', 'f8'), 
                           ('r', 'f8', 3), ('v', 'f8', 3)])
    
    # after changing the data type, which groups one row (8 numbers) to one type
    # we have one spare dimension (whose size is 1); let's squeeze it out
    data = np.squeeze(data)
    
    # return the required information
    return time, total, data