import numpy as np
from scipy.io import FortranFile
import quad_routs as qo

f = FortranFile('tunfw.dat','r')
i1 = f.read_ints(np.int32)
i2 = f.read_ints(np.int32)
i3 = f.read_ints(np.int32)
a = f.read_reals(float).reshape((3,4), order="F")
b = f.read_reals(float).reshape((5,2), order="F")
f.close()

iref = 2
iasp = 1
norder = 3
fname = './data/stell_norder'+str(norder)+'_iref'+str(iref)+'_iasp'+str(iasp)+'.bin'

s,t,nf = qo.read_surface_data(fname)
