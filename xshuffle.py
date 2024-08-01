#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Diffraction pattern resampling from stream for the analysis of coordinate errors
in TR-SFX crystallography data.

Requires
--------
    process.sh
    [OPTIONAL] parameters.yml

Usage
-----
1. Edit the variablesin the 'USER INPUT PARAMETERS' section as described below.
2. Run from the terminal:
    ./xshuffle.py --params parameters.yaml

References
----------
- Vallejos, A., Katona, G., Neutze, R.,
    Appraising protein conformational changes by resampling time-resolved 
    serial x-ray crystallography data
    Struct. Dyn. 11, 044302 (2024) https://doi.org/10.1063/4.0000258
"""

__version__     = "240724"
__author__      = "Adams Vallejos"
__email__       = "adams.vallejos.donoso <at> gu.se"
__maintainer__  = "Adams Vallejos"
__copyright__   = "Copyright 2024, Neutze Lab"
__license__     = "MIT license"
__status__      = "Production"
__credits__     = ["Adams Vallejos", "Gergely Katona","Richard Neutze"]


from re import search
from os import makedirs
from os.path import exists, split, splitext
from subprocess import check_call, call
from random import sample, choices
from time import localtime, strftime
import argparse
from numpy import ceil, int16
import yaml
from pandas import DataFrame

W = '\033[0m'
R = '\033[31m'
G = '\033[32m'
Y = '\033[33m'

# PDB PARSER
class Atom:
    '''
    Represents the 3D coordinates of a PDB 'ATOM' entry.
    
    '''
    def __init__(self, pdb_line):
        self.x_coordinate = float(pdb_line[30:38])
        self.y_coordinate = float(pdb_line[38:46])
        self.z_coordinate = float(pdb_line[46:54])

class PDB:
    '''
    Parses a PDB file to extract atom coordinates

    '''
    def __init__(self, file):
        self.file = file
        self.atoms = []
        self.parse()

    def parse(self):
        '''
        Parses a PDB file and extracts ATOM and HETATM records into Atom objects.

        '''
        f = open(self.file, 'r')
        for line in f.readlines():
            if line.startswith('ATOM') or line.startswith('HETATM'):
                atom = Atom(line)
                self.atoms.append(atom)
        f.close()

    def get_atoms(self):
        f = DataFrame.from_dict([x.__dict__ for x in self.atoms])
        return f


def read_arguments():
    '''
    Prints relevant information such as usage, version and more.

    '''
    parser = argparse.ArgumentParser(
        prog='xshuffle.py',
        description='''IMAGE RESAMPLING FROM STREAM FOR STRUCTURAL ANALYSIS OF TR-SX DATA.''',
        epilog=f'''{Y}"Draupnir, the dripper, was a gold ring made by the dwarves for Odin. Every
        ninth night, eight new rings of equal size and weight drip from it, but
        without the same magic ability as the original one."{W}'''
    )

    parser.add_argument('-V','--version',action='version', version=f' Version {__version__}')

    # mandatory = parser.add_argument_group(title='required arguments')
    # mandatory.add_argument("--params",required=True,
    #                        help='file containing custom processing parameters')
    parser.add_argument("--params",required=True,
                        help='file containing custom processing parameters')
    parser.add_argument("--coordinates",required=False, type=str, help='model file for phases')
    parser.add_argument("--stream",required=False, type=str, help='stream file to resample from')
    parser.add_argument("--fixed_size", required=False, type=int, help='Number of frames to select')
    parser.add_argument("--output", required=False, type=str, help='Path to store output files')

    return vars(parser.parse_args())


def timestamp():
    """
    Generate a string with current date and time

    """
    return strftime("%y-%m-%d %H:%M", localtime())

def welcome():
    ''' Prints a custom welcome message. '''
    print(f'''
    ############################################################
    ###  XSHUFFLE                             version {__version__}  ##
    ############################################################
    ''')


def verify_resampling(resampling_method):
    """
    Verify if the correct resampling resampling_method has been selected

    Parameters
    ==========
    resampling_method: str
        Either 'J' for Jackknife or 'B' for Bootstrap are allowed
    RESAMPLING_THRESH:  float
        Percentage of crystals to be used by Jackknife, default is None for
        Bootstrap.

    """
    if resampling_method[0].upper() == 'B':
        return print(timestamp(),f'Preparing {Y}Bootstrap{W} resampling...', flush=True)
    if resampling_method[0].upper() == 'J':
        return print(timestamp(),f'Preparing {Y}Jackknife{W} resampling...', flush=True)
    # if resampling_method[0].upper() == 'H':
    #     return print(timestamp(),f'Preparing {Y}50/50{W} resampling...', flush=True)

    raise TypeError(f'{R}Unrecognized resampling method {resampling_method} {W}')


def sample_from(chunks, resampling_method='B', resampling_thresh=None, fixed_size=None):
    '''
    Returns a dictionary with image number and number of copies used according
    to the resampling RESAMPLING_METHOD selected.

    Parameters
    ----------
    chunks: dict
        Images to resample
    resampling_method: str
        default value 'B' for bootstrap, also 'J' for Jackknife is available
    RESAMPLING_THRESH: float
        Default 'None' for bootstrap, if 'J' resampling_method is used, then give a
        RESAMPLING_THRESH value

    '''
    if resampling_method == 'J':
        if isinstance(fixed_size,int):
            k_len = fixed_size
        else:
            k_len = int(ceil(len(chunks)*resampling_thresh))
        j_sample = sorted(sample(list(chunks), k=k_len))
        j_dict = {i:j_sample.count(i) for i in j_sample}
        return j_dict
    if resampling_method == 'B':
        if isinstance(fixed_size,int):
            k_len = fixed_size
        else:
            k_len = len(chunks)
        b_sample = sorted(choices(list(chunks), k = k_len))
        b_dict = {i:b_sample.count(i) for i in b_sample}
        return b_dict
    if resampling_method == 'H':
        k_len = len(chunks)//2
        h_sample_a = sorted(sample(list(chunks), k=k_len))
        h_sample_b = [i for i in list(chunks) if i not in h_sample_a]

        h_dict_a = {i:h_sample_a.count(i) for i in h_sample_a}
        h_dict_b = {i:h_sample_b.count(i) for i in h_sample_b}
        return h_dict_a, h_dict_b

    raise ValueError(f'{R}Unsupported resampling method: {resampling_method}{W}')

def pdb_size(file, margin,spacing):
    atoms = PDB(file).get_atoms()
    limits = ceil([min(atoms.x_coordinate) - margin, max(atoms.x_coordinate) + margin,
                   min(atoms.y_coordinate) - margin, max(atoms.y_coordinate) + margin,
                   min(atoms.z_coordinate) - margin, max(atoms.z_coordinate) + margin])
    cellwork = int16([limits[1]-limits[0], limits[3]-limits[2], limits[5]-limits[4], 90, 90, 90])
    gridwork = int16(cellwork[:3]/spacing)
    xyzlim = int16(limits/spacing)

    return ' '.join(map(str, cellwork)), ' '.join(map(str, gridwork)),' '.join(map(str, xyzlim))

if __name__ == '__main__':

    welcome()

    arguments = read_arguments()

    # --- Read parameters from 'parameters.yaml'
    with open(arguments['params'],'r', encoding="utf-8") as p:
        params = yaml.load(p, Loader=yaml.FullLoader)
        ########################
        try:
            stream = params['stream']
        except KeyError:
            stream = None

        if arguments['stream']:
            stream = arguments['stream']
        if stream is None:
            raise TypeError(f'{R}Stream file missing{W}')
        ##########################
        try:
            coordinates = params['coordinates']
        except KeyError:
            coordinates = None

        if arguments['coordinates']:
            coordinates = arguments['coordinates']
        if coordinates is None:
            raise TypeError(f'{R}coordinates file missing{W}')
################################
        try:
            INPUT_DIR=params['input']
        except KeyError:
            INPUT_DIR = ''


        output = params['output']
        if arguments['output']:
            output = arguments['output']

        cell = params['cell']
        space_group = params['space_group']
        point_group = params['point_group']
        # coordinates = params['coordinates']
        library = params['cif_library']
        samples = params['samples']
        start = params['sample_start']
        RESAMPLING_METHOD = params['method'][0].upper()
        try:
            RESAMPLING_THRESH = params['threshold']
        except KeyError:
            RESAMPLING_THRESH = None

        try:
            FIXED_SIZE = params['fixed_size']
        except KeyError:
            FIXED_SIZE = None
        if arguments['fixed_size']:
            FIXED_SIZE = arguments['fixed_size']

        process_script = params['process_script']

    verify_resampling(RESAMPLING_METHOD)

    print(timestamp(), f'Reading file: {split(stream)[-1]}', flush=True)


    file_name = splitext(split(stream)[-1])[0]

    #   --- Stream cleaning ---
    event = {}
    TOTAL_IMAGES = 0
    INDEX_IMAGES = 0
    # with open(input_dir + stream, 'r', encoding="utf-8") as df:
    with open(stream, 'r', encoding="utf-8") as df:
        for line_number, line in enumerate(df):
            if line == '----- End chunk -----\n':
                event[line_number] = TOTAL_IMAGES
                TOTAL_IMAGES += 1
            elif line == '--- Begin crystal\n':
                INDEX_IMAGES += 1

    print(timestamp(),
        f'  {INDEX_IMAGES}/{TOTAL_IMAGES} chunks contain crystals.', flush=True)

    # Create 'output/'' folder if doesn't exist
    # print(output)
    if not exists(output):
        makedirs(output)

    print(timestamp(), 'Copying crystals from stream...', flush=True)
    COUNTER = 0
    INDEXED = False
    NUM_CHUNKS = 0
    NUM_INDEXED = 0
    crystals = {}
    # with open(input_dir + stream, 'r', encoding="utf-8") as data_1:
    with open(stream, 'r', encoding="utf-8") as data_1:
        # with open(input_dir + stream, 'r', encoding="utf-8") as data_2:
        with open(stream, 'r', encoding="utf-8") as data_2:

            # Create 'output/FILENAME' folder if doesn't exist
            streamdir = f'{output}{file_name}'
            if not exists(streamdir):
                makedirs(streamdir)

            # Reassign parameters after cleaning
            new_data = f'{streamdir}/clean.stream'

            with open(new_data,'w', encoding="utf-8") as clean_file:
                for line_number, line in enumerate(data_1):

                    if line_number == 0:
                        stream_header = line
                        clean_file.write(stream_header)

                    if search('indexed_by', line):
                        if search('none', line):
                            INDEXED = False
                        else:
                            INDEXED = True

                    if line == '----- End chunk -----\n':
                        NUM_CHUNKS += 1

                        if INDEXED:
                            NUM_INDEXED += 1
                            crystals[NUM_INDEXED] = line_number

                        while COUNTER <= line_number :
                            out = data_2.readline()
                            if INDEXED :
                                clean_file.write(out)

                            COUNTER += 1
                        INDEXED = False

        with open(f'{streamdir}/frame_loc.yaml','w', encoding="utf-8") as tag:
            yaml.dump(crystals, tag)

        print(timestamp(),
            f'  Indexed patterns copied to: {file_name}_clean.stream', flush=True)

    # --- CALCULATE PDB SIZE FOR MAPROT ---
    CELLWORK, GRIDWORK, XYZLIM = pdb_size(file=INPUT_DIR+coordinates,
                                          margin=5,
                                          spacing=0.25)

        # --- START RESAMPLING ---
    print(timestamp(),'Generating samples...', flush=True)


    # Select frames
    for current_sample in range(start, start + samples):
        sampledir = f'{streamdir}/{RESAMPLING_METHOD}'
        if not exists(sampledir):
            makedirs(sampledir)

        selected_chunks = sample_from(crystals,
                                      RESAMPLING_METHOD,
                                      RESAMPLING_THRESH,
                                      FIXED_SIZE)

        sample_path_str = f'{sampledir}/{current_sample:03}'

        if not exists(sample_path_str):
            makedirs(sample_path_str)

        # with open(f'{sample_path_str}/{current_sample:03}.yaml',
        with open(f'{sample_path_str}/nframes.yaml',
                  'w', encoding="utf-8") as sel:
            yaml.dump(selected_chunks, sel)

        with open(new_data, 'r', encoding="utf-8") as data_1:
            with open(f'{sample_path_str}/tmp.stream',
                      'w', encoding="utf-8") as sample_file:
                WRITE = False
                CHUNK = 0
                TEMPSEL = 0
                temp = []
                sample_file.write(stream_header)

                # print(sample_file.name)
                for line in data_1:
                    if line == '----- Begin chunk -----\n':
                        CHUNK += 1
                        if CHUNK in selected_chunks:
                            TEMPSEL += 1
                            WRITE = True

                    if WRITE :
                        temp.append(line)

                    if WRITE  and line == '----- End chunk -----\n':
                        WRITE = False

                        for i in range(selected_chunks[CHUNK]):
                            for item in temp:
                                sample_file.write(item)
                        temp = []

        print(timestamp(),
        f'{Y}-> Sample #{current_sample} generated!, processing{W}', flush=True)
        # Call 'process_script.sh' for processing
        check_call([f'./{process_script}',
        str(cell), # {0} $1
        str(space_group), # {1} $2
        str(point_group), # {2} $3
        str(sample_file.name), # {3} $4
        str(INPUT_DIR + coordinates), # {4} $5
        str(INPUT_DIR + library), # {5} $6
        str(CELLWORK), # {6} $7
        str(GRIDWORK), # {7} $8
        str(XYZLIM) # {8} $9
        ])

    # Run structural analysis and write out modified pdb

    # Draw the figure using 'resampling_render.pml'

    print(timestamp(), f'{G}TASK FINISHED!{W}', flush=True)
    call(['paplay','/usr/share/sounds/freedesktop/stereo/complete.oga'])
