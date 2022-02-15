import numpy as np
from scipy.io import FortranFile

class surface_info():
    sid = 0
    npatches = 0
    npts = 0
    srcvals = None
    srccoefs = None
    iptype = None
    norders = None
    ixyzs = None
    wts = None

class target_info():
    tid = 0
    ndtarg = 3
    ntarg = 0
    targets = None
    ipatch_id = None
    uvs_targ = None


class nearfield_info():
    rfac = 1.25
    nnz = 0
    nquad = 0
    sid = 0
    tid = 0
    row_ptr = None
    col_ind = None
    col_ptr = None
    row_ind = None
    iper = None
    iquad = None

def read_surface_data(fname):
    s = surface_info()
    t = target_info()
    nf = nearfield_info()
    s.sid = 0
    t.tid = 0
    nf.sid = s.sid
    nf.tid = t.tid

    f = FortranFile(fname,'r')
    s.npatches = int(f.read_ints(np.int32))
    s.npts = int(f.read_ints(np.int32))
    nf.nnz = int(f.read_ints(np.int32))
    t.ndtarg = int(f.read_ints(np.int32))
    t.ntarg = int(f.read_ints(np.int32))
    nf.nquad = f.read_ints(np.int32)
    s.norders = f.read_ints(np.int32)
    s.ixyzs = f.read_ints(np.int32)
    s.iptype = f.read_ints(np.int32)
    s.srcvals = f.read_reals(float).reshape((12,s.npts),order="F")
    s.srccoefs = f.read_reals(float).reshape((9,s.npts),order="F")
    s.wts = f.read_ints(np.int32)
    
    t.targets = f.read_reals(float).reshape((t.ndtarg,t.ntarg),order="F")
    t.ipatch_id = f.read_ints(np.int32)
    t.uvs_targ = f.read_reals(float).reshape((2,t.ntarg),order="F")
    
    nf.row_ptr = f.read_ints(np.int32)
    nf.col_ind = f.read_ints(np.int32)
    nf.iquad = f.read_ints(np.int32)
    nf.col_ptr = f.read_ints(np.int32)
    nf.row_ind = f.read_ints(np.int32)
    nf.iper = f.read_ints(np.int32)
    f.close()
    return s,t,nf
    
    
    



