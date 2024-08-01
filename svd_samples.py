#!/usr/bin/env python3
###############################################################
### svd_samples                             version 240724   ##
### (C) 2024 Adams Vallejos                                  ##
###############################################################

import shutil
import mrcfile
import sys
import numpy
from numpy import pad, zeros, sqrt, matmul, matrix, array
from numpy.linalg import eig, svd
import os
from os.path import basename, splitext
from pandas import DataFrame, set_option
import matplotlib.pyplot as plt

# Use a class to retrieve all information in a table, the labels are the following
class MrcHeader:
    def __init__(self):
        self.nx = None
        self.ny = None
        self.nz = None
        self.mode = None
        self.nxstart = None
        self.nystart = None
        self.nzstart = None
        self.mx = None
        self.my = None
        self.mz = None
        self.cella = None
        self.cellb = None
        self.mapc = None
        self.mapr = None
        self.maps = None
        self.dmin = None
        self.dmax = None
        self.dmean = None
        self.ispg = None
        self.nsymbt = None
        self.extra1 = None
        self.exttyp = None
        self.nversion = None
        self.extra2 = None
        self.origin = None
        self.map = None
        self.machst = None
        self.rms = None
        self.nlabl = None
        self.label = None
        self.data = None
        self.voxel_size = None
        
    def parse(self, filename):
        with mrcfile.open(filename,'r') as mrc:
            # TODO: Set the proper type to every attribute
            self.nx = int(mrc.header.nx)
            self.ny = int(mrc.header.ny)
            self.nz = int(mrc.header.nz)
            self.mode = int(mrc.header.mode)
            self.nxstart = int(mrc.header.nxstart)
            self.nystart = int(mrc.header.nystart)
            self.nzstart = int(mrc.header.nzstart)
            self.mx = int(mrc.header.mx)
            self.my = int(mrc.header.my)
            self.mz = int(mrc.header.mz)
            self.cella = array(mrc.header.cella)
            self.cellb = mrc.header.cellb
            self.mapc = int(mrc.header.mapc)
            self.mapr = int(mrc.header.mapr)
            self.maps = int(mrc.header.maps)
            self.dmin = mrc.header.dmin
            self.dmax = mrc.header.dmax
            self.dmean = mrc.header.dmean
            self.ispg = mrc.header.ispg
            self.nsymbt = mrc.header.nsymbt
            self.extra1 = mrc.header.extra1
            self.exttyp = mrc.header.exttyp
            self.nversion = mrc.header.nversion
            self.extra2 = mrc.header.extra2
            self.origin = mrc.header.origin
            self.map = mrc.header.map
            self.machst = mrc.header.machst
            self.rms = mrc.header.rms
            self.nlabl = mrc.header.nlabl
            self.label = mrc.header.label
            self.data = mrc.data
            self.voxel_size = mrc.voxel_size
    
    def get_data(self):
        return self.data
        
    def to_dict(self):
        return {
            'nx': self.nx,
            'ny': self.ny,
            'nz': self.nz,
            'mode': self.mode,
            'nxstart': self.nxstart,
            'nystart': self.nystart,
            'nzstart': self.nzstart,
            'mx': self.mx,
            'my': self.my,
            'mz': self.mz,
            'cella': self.cella,
            'cellb': self.cellb,
            'mapc': self.mapc,
            'mapr': self.mapr,
            'maps': self.maps,
            'dmin': self.dmin,
            'dmax': self.dmax,
            'dmean': self.dmean,
            'ispg': self.ispg,
            'nsymbt': self.nsymbt,
            'extra1': self.extra1,
            'exttyp': self.exttyp,
            'nversion': self.nversion,
            'extra2': self.extra2,
            'origin': self.origin,
            'map': self.map,
            'machst': self.machst,
            'rms': self.rms,
            'nlabl': self.nlabl,
            'label': self.label,
            'voxel_size': self.voxel_size,
        }
# END OF HEADER CLASS

# Map dimensions nx, ny, nz These vary from map to map
# corresponding to map dimensions NC, NR, NS (fast axis, medium axis, slow axis)
# These are orientated According to MAPC, MAPR, MAPS keywords

# WHEN TAKING THE DATA FROM AN ORIGINAL MAP USE MODE 'r' ONLY

# START BY READING THE HEADERS
def set_workspace(workspace):
    if not os.getcwd() == workspace:
        os.chdir(workspace)
    return f'Working on: {workspace}'


def maps_from_list(list):
    files = []
    with open(list) as lst:
        for item in lst:
            files.append(item.strip('\n'))
    return files


def block_dimensions(list):
    '''
    Returns a dictionary of map of {axis: max dimension} from a list of ccp4 maps 
    
    list: lst
        list of files to read    
    
    Map dimensions NX, NY, NZ vary from map to map, 
    corresponding to map dimensions NC, NR, NS (fast axis, medium axis, slow axis).
    
    These are orientated According to MAPC, MAPR, MAPS keywords
    '''
    nc = 0
    nr = 0
    ns = 0
    for i in range(len(list)):
        # Close file after reading data to avoid overloading memory
        with mrcfile.open(list[i],mode='r',header_only=True) as f:
            nc = max(int((nc)), int(f.header.nx))
            nr = max(int((nr)), int(f.header.ny))
            ns = max(int((ns)), int(f.header.nz))
            # Switch dimension labels 1,2,3 -> 0,1,2
            mapc = f.header.mapc - 1
            mapr = f.header.mapr - 1
            maps = f.header.maps - 1
    return dict([(mapc, nc), (mapr, nr), (maps, ns)])

def check_origin(list):
    '''
    Verifies if the origin of the coordinate system of all maps is to (0, 0, 0).
    
    
     Example:
        rec.array((0., 0., 0.),dtype=[('x', '<f4'), ('y', '<f4'), ('z', '<f4')])
    
    Analogous to pandas table:
        - Keys: 'x', 'y', 'z'
        - type '<f4': Little endian 32 byte float, the least significant byte 
                      (the "little end") of the data is placed at the byte 
                      with the lowest address. 
    '''
    for i in range(len(list)):
        # Close file after reading data to avoid overloading memory
        with mrcfile.open(list[i],mode='r',header_only=True) as f:
            if i == 0:
                origin = f.header.origin
                
            if f.header.origin != origin:
                raise TypeError('Unequal origin')
            
    if origin.item() == (0,0,0):
        return origin
    raise TypeError('Origin not at (0, 0, 0)')
    


# Load data from list
def load_data(list, as_df = True):
    '''
    Returns 
        data: A list of 3D arrays from a list of ccp4 maps 
        headers: A pandas dataframe with information of all headers 
    
    '''
    data = []
    headers = []
    for mrc_file in list:
        header = MrcHeader()
        header.parse(mrc_file)
        
        headers.append(header.to_dict())
        data.append(header.get_data())
    
    if as_df:
        headers = DataFrame(headers)
    
    return headers, data



# Generate global variables according to header info:
def get_grid_specs(headers,swap_axes=None):
    '''
    Returns a dictionary of 'axis': units 
    and sets a number of global variables to preserve across all maps.
    
    '''
    # set the units to be preserved along all the maps
    
    # Select the unit cell intervals from the first as a reference
    if not headers.mx.eq(headers.mx[0]).all():
        raise TypeError('Unequal sampling grid along X axis!')
    
    if not headers.my.eq(headers.my[0]).all():
        raise TypeError('Unequal sampling grid along Y axis!')
    
    if not headers.mz.eq(headers.mz[0]).all():
        raise TypeError('Unequal sampling grid along Z axis!')
            
    MX, MY, MZ = headers.iloc[0]['mx'], headers.iloc[0]['my'], headers.iloc[0]['mz']
    globals()['MX'], globals()['MY'], globals()['MZ'] = MX, MY, MZ
    
    # Set the largest array dimension
    NX, NY, NZ = headers['nx'].max(), headers['ny'].max(), headers['nz'].max()
    globals()['NX'], globals()['NY'], globals()['NZ'] = NX, NY, NZ
    
    print(f'These values will be output to the 3D grid:\n'
          + f'nx: {NX}, ny: {NY}, nz: {NZ}, \n'
          + f'mx: {MX}, my = {MY}, mz = {MZ} \n')
    
    if swap_axes:
        # This is added since (mapc, mapr, mapc) in maps exported 
        # from different programs do not do not give the proper 
        # ordering of the data array as described in the header.
        return {int(swap_axes[i]) - 1: dim for i, dim in enumerate([NX, NY, NZ])}
        
    return {headers.iloc[0]['mapc'] - 1: NX, 
            headers.iloc[0]['mapr'] - 1: NY, 
            headers.iloc[0]['maps'] - 1: NZ}


# PERFORM ARRAY PADDING
def pad_array(list_of_arrays):
    padded_data = []
    for arr in list_of_arrays:
        # PAD WITH ZEROS AT THE START OF ARRAY
        # pad_width = [(map_specs[i]-arr.shape[i],0) for i in range(3)]
        
        # PAD WITH ZEROS AT THE END OF ARRAY
        pad_width = [(0,map_specs[i]-arr.shape[i]) for i in range(3)]
        padded_data.append(pad(array=arr, pad_width=pad_width, mode='constant', constant_values=0))
    return padded_data


# WRITE THE DATA TO OUTPUT COPY OF REFERENCE
def write_map(new_data, in_map, out_map):
    shutil.copy(in_map, out_map)
    with mrcfile.open(out_map, 'r+') as f:
        f.set_data(new_data)
        f.update_header_from_data()
        f.update_header_stats()
        f.header.mx = MX
        f.header.my = MY
        f.header.mz = MZ


#   --- USER INPUT PARAMETERS

LIGHT = "light"

RES='B'

WORKSPACE = f'/PATH/TO/RESAMPLING/output/{LIGHT}/{RES}'

out_dir = 'svd_FoFo'

print(f'Working on {WORKSPACE}')

set_workspace(WORKSPACE)

# Load a list of map data
# files = maps_from_list('../maps_fofo.lst')
files = []
for i in range(0,100):
    files.append(f'{i:03}/analysis_fofo_0/merge_trunc_cad_uniq_fofo_map_coeffs.ccp4')

print(f'Reading maps grid specs')

try:
    os.mkdir(out_dir)
    print(f"Directory '{out_dir}' created ")
except FileExistsError:
    print(f"Directory '{out_dir}' already exists")  

headers, data = load_data(files)

# Set the values to be used across the maps 
map_specs = get_grid_specs(headers,swap_axes=[3,2,1])
# map_specs

# PERFORM PADDING ACCORDING TO BLOCK SPECS
padded_data = pad_array(data)

print(f'Padding data arrays')
# DO SVD

# Generate a matrix of zeros of dimensions MAP SIZE x NUMBER OF MAPS
A = zeros((padded_data[0].size, len(files)))

for i in range(len(files)):
    A[:,i] = padded_data[i].reshape(-1)


print(f'Calculating SVD')
# SVD method_1
cov = A.T @ A 
L,V = eig(cov)

S1 = sqrt(L[0])
for i in range(A.shape[1]):
    print(f'{i+1}/{A.shape[1]}', end='')
    U1 = matmul(A,V[:,i])/S1
    Ar = S1*(matrix(U1).T)*(matrix(V[:,i]))
#     shutil.copy(files[i], f'{outdir}/{splitext(files[i])[0]}.ccp4')    
    S1 = sqrt(L[0])
    U1 = matmul(A,V[:,0])/S1
    Ar = S1*(matrix(U1).T)*(matrix(V[:,0]))
    Ar = array(Ar).astype(numpy.float32)

    # Copy a reference map and set the new data to it.
    write_map(Ar[:,0].reshape(padded_data[0].shape), in_map=files[i], out_map=f'{out_dir}/svd_{splitext(basename(files[i]))[0]}.ccp4')
    print('', end='\r')

print('JOB FINISHED!')

