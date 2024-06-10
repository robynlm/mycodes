"""
Author: Robyn Munoz

"""

import h5py
import sys
import re
import numpy as np
from tools import ReadingTools as RRead

class MyClass:
    def __init__(self, param, f):
        self.param = param
        self.iG = self.param['ghost_size']
        self.f = f
        self.it_filename = param['h5datapath']+param['simname']
        
    def cp_it(self, key):
        # Extract the iteration
        it = int(re.split(' |=', key)[2])
        # Copy that key to the iteration file
        fname = '{}_it_{:06d}.hdf5'.format(self.it_filename, it)
        fnew = h5py.File(fname,'a')
        fnew_keys = [k for k in fnew.keys()]
        logical = [k in key for k in fnew_keys]
        if np.sum(logical)==0:
            self.f.copy(key,fnew)
            # TODO copy time attribute
            
    def join_chunks(self):
        for rl in range(param['max_refinement_levels']):
            file_keys = [k for k in self.f.keys() if 'rl='+str(rl) in k]
            new_file_keys = list(set([k.split(' c=')[0] 
                                      for k in file_keys]))
            for new_key in new_file_keys:
                
                # number of sections
                nbrmpi = 0
                ghosts_on = False
                for k in file_keys:
                    if (new_key in k) and (' c=' in k):
                        nbrmpi = np.max([nbrmpi, int(k.split(' c=')[1])])
                        ghosts_on = True
                        
                # merge the sections
                if nbrmpi>0:
                    nbrmpi += 1
                    
                    cut_keys = [new_key+' c='+str(c) 
                                for c in range(nbrmpi)]
                    if ghosts_on:
                        cut_data = [np.array(self.f[k])[self.iG:-self.iG,
                                                        self.iG:-self.iG,
                                                        self.iG:-self.iG] 
                                    for k in cut_keys]

                    # =================
                    if nbrmpi == 1:
                        uncut_data = cut_data[0]
                    elif nbrmpi == 2:
                        uncut_data = np.append(
                            cut_data[0], cut_data[1], axis=0)
                    elif nbrmpi == 3:
                        uncut_data = np.append(
                            np.append(cut_data[0], cut_data[1], axis=1), 
                            cut_data[2], axis=0)

                    # =================
                    elif nbrmpi >= 4 and nbrmpi < 9:
                        # fix axis = 2
                        if nbrmpi == 4:
                            ndata = cut_data
                        elif nbrmpi == 5:
                            ndata = [np.append(
                                cut_data[0], cut_data[1], axis=2),
                                     cut_data[2], cut_data[3], cut_data[4]]
                        elif nbrmpi == 6:
                            ndata = [np.append(
                                cut_data[0], cut_data[1], axis=2), 
                                     cut_data[2],
                                     np.append(
                                         cut_data[3], cut_data[4], axis=2), 
                                     cut_data[5]]
                        elif nbrmpi == 7:
                            ndata = [np.append(
                                cut_data[i], cut_data[i+1], axis=2) 
                                     for i in np.arange(3)*2]
                            ndata += [cut_data[6]]
                        else:
                            ndata = [np.append(
                                cut_data[i], cut_data[i+1], axis=2) 
                                     for i in np.arange(4)*2]
                        # fix axis = 1 and 0
                        uncut_data = np.append(
                            np.append(ndata[0], ndata[1], axis=1), 
                            np.append(ndata[2], ndata[3], axis=1), 
                            axis=0)

                    # =================
                    elif nbrmpi >= 9 and nbrmpi < 13:
                        # fix axis = 2
                        if nbrmpi == 9:
                            ndata = [
                                np.append(
                                    cut_data[0], cut_data[1], axis=2), 
                                cut_data[2], 
                                np.append(
                                    cut_data[3], cut_data[4], axis=2), 
                                cut_data[5], 
                                np.append(
                                    cut_data[6], cut_data[7], axis=2), 
                                cut_data[8]]
                        elif nbrmpi == 10:
                            ndata = [
                                np.append(
                                    cut_data[0], cut_data[1], axis=2), 
                                np.append(
                                    cut_data[2], cut_data[3], axis=2), 
                                np.append(
                                    cut_data[4], cut_data[5], axis=2), 
                                cut_data[6],
                                np.append(
                                    cut_data[7], cut_data[8], axis=2), 
                                cut_data[9]]
                        elif nbrmpi == 11:
                            ndata = [
                                np.append(
                                    cut_data[i], cut_data[i+1], axis=2) 
                                for i in np.arange(5)*2]
                            ndata += [cut_data[10]]
                        else:
                            ndata = [
                                np.append(
                                    cut_data[i], cut_data[i+1], axis=2) 
                                for i in np.arange(6)*2]
                        # fix axis = 1 and 0
                        nndata = [np.append(ndata[i], ndata[i+1], axis=1) 
                                  for i in np.arange(3)*2]
                        uncut_data = np.append(
                            np.append(nndata[0], nndata[1], axis=0), 
                            nndata[2], axis=0)

                    # =================
                    else:
                        print('idk ze mpi :(', flush=True)
                        break

                    try:
                        self.f.create_dataset(new_key, data=uncut_data)
                        for atk in ['delta', 'time', 'level', 'timestep']:
                            self.f[new_key].attrs[atk] = self.f[cut_keys[0]].attrs[atk]
                    except:
                        print("already exists", flush=True)

                    for k in cut_keys:
                        del self.f[k]
        
        

if __name__ == "__main__":
    print("\n#################################", flush=True)
    print("\n      Split & Merge Files        \n", flush=True)
    print("#################################\n", flush=True)
    
    simname = sys.argv[1]
    param = RRead.read_parameters(simname)
    RRead.MakeDir(param['h5datapath'][:-1])
    
    print("\n========================== Split the files per iteration", 
          flush=True)  
    
    for irestart in range(param['nbr_restarts']):
        if irestart>=1:
            print("\n=============", flush=True)  
            print('Going through restart nbr '+str(irestart), flush=True)
        datapath = (param['HorSpath'] + param['simname'] 
                    + '/output-{:04d}/'.format(irestart) 
                    + param['simname'] + '/')
        # Get .h5 files to split
        files = RRead.BASH('ls ' + datapath + '*.h5').split('\n')  
        # Files as all HDF5 output by Cactus (excludes chkpts)
        print('files to be split: ', flush=True)
        files_kept = []
        last_part = []
        for file in files:
            if ('checkpoint' not in file) and ('NaN' not in file):
                files_kept += [file]
                last_part += [file.split('/')[-1]]
        print(last_part, flush=True)

        # Splitting the files per iteration
        want_to_copy_things = True
        if want_to_copy_things:
            for file in files_kept:
                if 'ERROR' in file:
                    print('.h5 file not present in :' 
                          + datapath, flush=True)
                else:
                    # Open the file
                    f = h5py.File(file,'r')
                    print('\nusing file: {}'.format(f.filename), 
                          flush=True)

                    # Collect the keys
                    filekeys = [key for key in f.keys() if 'it=' in key]
                    print('filekeys obtained', flush=True)
                    print(filekeys[0], ', ', 
                          filekeys[1], ' ... ', 
                          filekeys[-2], ', ', 
                          filekeys[-1], flush=True)
                    
                    # Check if file is already copied
                    command = 'ls '+param['h5datapath']+'*.hdf5'
                    files_already_done = RRead.BASH(command)
                    # Have any been copied so far?
                    if 'ERROR' in files_already_done:
                        copy_file = True
                    else:
                        last_it = np.max([int(re.split(' |=', k)[2]) 
                                          for k in filekeys])
                        last_file_name = (param['h5datapath'] 
                                          + param['simname'] 
                                          + '_it_' 
                                          + '{:06d}.hdf5'.format(last_it))
                        # is the last key in the last file?
                        if last_file_name in files_already_done:
                            last_file = h5py.File(last_file_name, 'r')
                            last_file_keys = [k for k in last_file.keys()]
                            last_keys = [k for k in filekeys 
                                        if 'it='+str(last_it) in k]
                            logical = [k in last_file_keys 
                                       for k in last_keys]
                            # have these keys already been done?
                            if np.sum(logical)==len(logical):
                                copy_file = False
                                print('already copied', flush=True)
                            else:
                                copy_file = True
                            last_file.close()
                            del last_keys, last_file
                            del last_file_keys, logical
                        else:
                            copy_file = True
                        del last_it, last_file_name
                    
                    copy_file = True
                    # Copy and split data
                    nbr_filekeys = len(filekeys)
                    if copy_file:
                        C = MyClass(param, f)
                        for keycounter, key in enumerate(filekeys):
                            C.cp_it(key)
                            if keycounter % int(nbr_filekeys/10) == 0:
                                percent_done = (keycounter * 100 
                                                / nbr_filekeys)
                                if keycounter==0:
                                    print('Progress = '
                                          + '{:.2f}%'.format(percent_done), 
                                          end="", flush=True)
                                else:
                                    print(', {:.2f}%'.format(percent_done), 
                                          end="", flush=True)
                                    
                        del C
                    print(flush=True)
                    f.close()
                    #RRead.BASH('rm '+file)
        
    # Merge mpi sections
    print("\n========================== Merge domain separated by mpi procs", 
          flush=True)
    
    command = 'ls ' + param['h5datapath'] + param['simname'] + '*'
    all_it_file_names = RRead.BASH(command).split('\n')
    nbr_it = len(all_it_file_names)
    
    for i, file in enumerate(all_it_file_names):
        f = h5py.File(file, "a")
        C = MyClass(param, f)
        C.join_chunks()
        f.close()
        del C
        percent_done = (i * 100 / nbr_it)
        print(re.split('iterations/|.hdf5', file)[1]
              + ' Progress = {:.2f}%'.format(percent_done), flush=True)
        
    
    print('Done!', flush=True)
