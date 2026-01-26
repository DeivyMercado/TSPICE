"""
tSPICE - Tidal Signal with Python and SPICE
A package for calculating the tidal potential using SPICE.
"""

#Version of the TSPICE package
__version__ = "0.0.2"

import os
import json
import spiceypy as spy

import urllib.request
import time

# To track if kernels are loaded
ROOTDIR = os.path.dirname(os.path.abspath(__file__))
kernels_loaded = False

def get_data(filename):
    """
    Get the full path of the `filename` which is one of the datafiles provided with the package.
    
    Parameters
    ----------
    filename : str
        Name of the data file.
        
    Returns
    -------
    path : str
        Full path to package datafile.
    """
    return os.path.join(ROOTDIR, 'data', filename)

def initialize(data_directory=None, verbose=True, level='minimal'):
    """
    Download and load the necessary SPICE kernels.
    
    Input:
    - data_directory : [str, optional] Path to the directory where kernels should be stored and loaded from. 
                       If not provided, a 'tspice_data' directory will be created in the current working directory.
    - verbose : [bool, optional] If True, print detailed information about kernel downloading and loading. Default is True.
    - level : [str, optional] Level of kernels to download. Options are 'minimal', 'planetary', 'full'. Default is 'full'.
              'minimal': Only small kernels (non-.bsp).
              'planetary': Minimal + Planetary ephemeris (de442.bsp).
              'full': All kernels.
    """
    global kernels_loaded
    
    # Define data directory
    if data_directory is None:
        data_dir = os.path.join(os.getcwd(), 'tspice_data')
    else:
        data_dir = data_directory

    # --- logic from download_kernels ---
    #Directory for the kernels and the configuration file
    # kernel_dir = os.path.join(data_dir, 'kernels')
    # User requested to avoid using 'kernels' subdirectory, so using data_dir directly
    kernel_dir = data_dir 
    
    # kernel_config is still inside the package installation
    kernel_config = os.path.join(ROOTDIR, 'data', 'spice_kernels.json')

    #Read the JSON file
    with open(kernel_config, 'r') as f:
        config = json.load(f) 	#Reads the file and turns it into a Python dictionary!
    
    # Filter kernels based on level
    levels_map = {
        'minimal': ['minimal'],
        'planetary': ['minimal', 'planetary'],
        'full': ['minimal', 'planetary', 'full']
    }
    
    if level not in levels_map:
        raise ValueError(f"Invalid level '{level}'. Must be one of {list(levels_map.keys())}.")
        
    allowed_levels = levels_map[level]
    kernels_to_download = [k for k in config['kernels'] if k.get('level', 'full') in allowed_levels]

    #If the kernel directory does not exist, create it and download the kernels
    #Note: initialize creates data_dir implicitly if we pass it to makedirs, or valid checks
    if not os.path.exists(kernel_dir):
        os.makedirs(kernel_dir)
        print(f"Directory for Kernels created at: {kernel_dir}")
        
        num_files = len(kernels_to_download)
        print(f"WARNING: You are about to download {num_files} kernel files (level='{level}').")
        print("This process may take a long time depending on your internet connection.")
        
        #Start downloading the kernel files
        for i,k in enumerate(kernels_to_download):
            filename = k['filename']
            url = k['url']
            filepath = os.path.join(kernel_dir, filename)
            part_filepath = filepath + '.part'
            
            print(f"\nDownloading file {i+1} of {num_files}: {filename}")
            if verbose: print(f"Source: {url}")
            
            try:
                # Progress bar reporthook
                def reporthook(count, block_size, total_size):
                    if total_size > 0:
                        percent = int(count * block_size * 100 / total_size)
                        # Simple text progress bar
                        width = 50
                        filled = int(width * percent / 100)
                        bar = '=' * filled + '-' * (width - filled)
                        print(f"\r[{bar}] {percent}%", end='')
                
                urllib.request.urlretrieve(url, part_filepath, reporthook=reporthook)
                print() # Newline after progress bar finishes
                
                # Rename .part to actual filename upon success
                os.rename(part_filepath, filepath)
                if verbose: print(f"Successfully downloaded {filename}.")
                
            except Exception as e:
                print(f"\nError downloading the file {filename} from {url}.")
                print(f"Exception: {e}")
                # Optional: clean up partial file
                if os.path.exists(part_filepath):
                    os.remove(part_filepath)

    #If the directory exists, verify that the files are there already
    else:
        if verbose: print(f"Directory for Kernels already exists at: {kernel_dir}")
        
        #Check if we need to verify/resume downloads? For now, we stick to the basic check
        #Verify that the files are there already or download them if they are not
        missing_kernels = []
        for k in kernels_to_download:
            if k['filename'] not in os.listdir(kernel_dir):
                missing_kernels.append(k)
        
        if missing_kernels:
            print(f"Found {len(missing_kernels)} missing kernels for level '{level}'. improved download starting...")
            num_files = len(missing_kernels)
            
            for i,k in enumerate(missing_kernels):
                filename = k['filename']
                url = k['url']
                filepath = os.path.join(kernel_dir, filename)
                part_filepath = filepath + '.part'
                
                print(f"\nDownloading missing file {i+1} of {num_files}: {filename}")
                
                try:
                     # Progress bar reporthook
                    def reporthook(count, block_size, total_size):
                        if total_size > 0:
                            percent = int(count * block_size * 100 / total_size)
                            width = 50
                            filled = int(width * percent / 100)
                            bar = '=' * filled + '-' * (width - filled)
                            print(f"\r[{bar}] {percent}%", end='')
                    
                    urllib.request.urlretrieve(url, part_filepath, reporthook=reporthook)
                    print()
                    
                    # Rename .part to actual filename upon success
                    os.rename(part_filepath, filepath)
                    if verbose: print(f"Successfully downloaded {filename}.")

                except Exception as e:
                    print(f"\nError downloading {filename}: {e}")
                    if os.path.exists(part_filepath):
                        os.remove(part_filepath)
        else:
             if verbose: print(f"All kernels for level '{level}' appear to be present.")

    
    # --- logic from write_meta_kernel ---
    #Create the path to the meta kernel
    meta_kernel_path = os.path.join(data_dir, 'meta_kernel')

    #If the meta kernel does not exist, create it, or update it if level changed
    # We simply overwrite it every time to ensure it matches the requested level
    with open(meta_kernel_path, 'w') as f:

        f.write('KPL/MK\n')	#Type of file
        f.write('\n')	#Beginning of data
        f.write('\\begindata\n')	#Beginning of data

        f.write(f"PATH_VALUES = ('{kernel_dir}')\n")	#Path to the kernels
        f.write("PATH_SYMBOLS = ('KERNELS')\n")	#Symbol to indicate the path
        f.write('\n')	#New line	

        #Writing the kernels
        f.write('KERNELS_TO_LOAD = (\n')
        for k in kernels_to_download:
            filename = k['filename']
            f.write(f"'$KERNELS/{filename}',\n")
        f.write(')\n')
    
    if verbose: print(f"Meta kernel created/updated at {meta_kernel_path}.")
    
    #Load the meta kernel
    spy.furnsh(meta_kernel_path)
    
    kernels_loaded = True
    print(f"TSPICE initialized successfully. Kernels loaded from: {meta_kernel_path}")
    if verbose: print(f"Kernels loaded from: {meta_kernel_path}")


#Import key functions to make them available at the top level for users. This allows users to write:
#from tspice import initialize, Body, BodyResponse
#Instead of:
#from tspice.tidal_signal import initialize, Body
#from tspice.internal_dynamics import BodyResponse

from .tidal_signal import Body
from .internal_dynamics import BodyResponse
from . import utils
from . import integration_tools
from .plotting import plot_one_signal, plot_many_signal
# from .data import Data


# Package initialization message
def _welcome_message():
    """Display welcome message on import."""
    print(f"Welcome to tSPICE v{__version__}")

_welcome_message()

# Clean up namespace
del _welcome_message

#List what should be imported with "from tspice import *"
__all__ = ['initialize', 'Body', 'BodyResponse', 'internal_dynamics', 'utils', 'integration_tools', 'plot_one_signal', 'plot_many_signal', 'get_data'] #, 'Data']