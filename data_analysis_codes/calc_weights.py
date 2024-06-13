import os
import h5py
import sys
import numpy as np
from tools.lineint import lineint
from tools import ReadingTools as RRead
from multiprocessing import Pool, Lock, Value

#lock = Lock()
# only one process writes to the file at a time
"""
def init(args):
    global it_counter
    it_counter = args
"""

def calc_things(li, f, ix, iy, iz):
    end_point = (ix, iy, iz)
    outputs = li.segments_idx_end(end_point)
    indices, segments, theta, phi = outputs
    idx = '({},{},{})'.format(ix, iy, iz)
    #with lock:
    f.create_dataset(idx+' indices', data=indices)
    f.create_dataset(idx+' segments', data=segments)
    """
    # Progress
    global it_counter
    with it_counter.get_lock():
        it_counter.value += 1
    percent_done = (it_counter.value * 100 / nbr_iterations)
    print('Progress = {:.2f}%'.format(percent_done), flush=True)
    """
        
if __name__ == "__main__":
    
    simname = sys.argv[1]
    #nbr_proc = int(sys.argv[2])
    
    param = RRead.read_parameters(simname)
    print('', flush=True)
    print('Simulation:', param['simname'], flush=True)
    rl = param['max_refinement_levels'] - 1
    
    # size smallest resolution grid
    it = 0
    f = h5py.File(param['h5datapath'] + param['simname'] 
                  + '_it_{:06d}.hdf5'.format(it),'r')
    key = 'ADMBASE::alp it={} tl=0 rl={}'.format(it, rl)
    if rl > 1:
        Nsg = np.shape(np.array(f[key])[3:-3, 3:-3, 3:-3])[0]
    else:
        Nsg = np.shape(np.array(f[key]))[0]
    f.close()
    print('number of grid points on smallest resolution grid:', Nsg, flush=True)
    
    # def param of smallest resolution grid
    dxg = param['dx'] / (2**rl)
    Lg = Nsg * dxg
    xmin = (- Lg / 2 + (dxg / 2) * (Nsg % 2))
    
    # define lineint
    xdomain = [xmin, xmin + param['Lx'], dxg]
    Ng = param['Lx'] / dxg
    starting_point = [int(Ng/2)] * 3
    li = lineint.LineIntegrate(xdomain, xdomain, xdomain, starting_point)
    print('lineint defined', flush=True)
    
        
    # remove file if exists
    filename = '/users/munozr/simulations/' + param['simname'] + '/proper_distance_weight.hdf5'
    if os.path.exists(filename):
        os.remove(filename)
        print(filename, 'removed', flush=True)

    # calc weight and save
    """
    print(' === ready to rumble', flush=True)
    f = h5py.File(filename, 'w') 
    with Pool(processes=nbr_proc) as pool:
        for ix in ixall:
            for iy in ixall:
                for iz in ixall:
                    pool.apply_async(calc_things, 
                                     args=(li, f, ix, iy, iz, nbr_iterations))
        pool.close()
        pool.join()
    """
    print(' === ready to rumble', flush=True)
    f = h5py.File(filename, 'w') 
    it = 0
    nbr_iterations = li.Nx*li.Ny
    iz = int(li.Nz/2)
    for ix in range(li.Nx):
        for iy in range(li.Ny):
            calc_things(li, f, ix, iy, iz)
            it += 1
        percent_done = (it * 100 / nbr_iterations)
        print('Progress = {:.2f}%'.format(percent_done), flush=True)
    f.close()