##################################################################
#                                                                #
# tSPICE: Tidal Signal with Python and SPICE                    #
#                                                                #
##################################################################
# License: GNU Affero General Public License v3 (AGPL-3.0)       #
# Authors: Jorge I. Zuluaga, Juanita A. Agudelo                  #
# Contact: jorge.zuluaga@udea.edu.co                             #
##################################################################
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
    Get the full path to a package data file.
    
    This function returns the absolute path to data files included with the tSPICE package,
    such as SPICE kernel configuration files or planetary models.
    
    Parameters
    ----------
    filename : str
        Name of the data file (e.g., 'spice_kernels.json', 'PREM_amorin.xlsx').
        
    Returns
    -------
    str
        Absolute path to the requested data file.
        
    Examples
    --------
    >>> import tspice
    >>> kernel_config = tspice.get_data('spice_kernels.json')
    >>> print(kernel_config)
    /path/to/tspice/data/spice_kernels.json
    """
    return os.path.join(ROOTDIR, 'data', filename)

def initialize(data_directory=None, verbose=True, level='planetary'):
    """
    Download and load SPICE kernels for tSPICE.
    
    This function handles downloading SPICE kernels from NASA's NAIF server and loading them
    into memory. Kernels are essential for computing ephemerides of celestial bodies. The function
    creates a meta-kernel file that references all downloaded kernels and loads them using SPICE.
    
    Parameters
    ----------
    data_directory : str, optional
        Path to the directory where kernels should be stored and loaded from.
        If not provided, a 'tspice_data' directory will be created in the current working directory.
    verbose : bool, optional
        If True, print detailed information about kernel downloading and loading.
        Default is True.
    level : {'minimal', 'planetary', 'full'}, optional
        Level of kernels to download. Default is 'planetary'.
        
        - 'minimal': Only small kernels (leap seconds, frame definitions, etc.)
        - 'planetary': Minimal + planetary ephemeris (de442.bsp, ~600 MB)
        - 'full': All kernels including satellite ephemerides
    
    Raises
    ------
    ValueError
        If an invalid `level` argument is provided.
    
    Notes
    -----
    Downloaded kernels are cached in the specified directory. If kernels already exist,
    they will not be re-downloaded unless they are missing.
    
    The function sets the global `kernels_loaded` flag to True upon successful completion.
    
    Examples
    --------
    Initialize with default settings (planetary kernels in current directory):
    
    >>> import tspice
    >>> tspice.initialize()
    
    Initialize with minimal kernels in a custom directory:
    
    >>> tspice.initialize(data_directory='/path/to/kernels', level='minimal')
    
    Quiet initialization:
    
    >>> tspice.initialize(verbose=False)
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
    

    spy.furnsh(meta_kernel_path)
    
    kernels_loaded = True
    print(f"TSPICE initialized successfully. Kernels loaded from: {meta_kernel_path}")
    if verbose: print(f"Kernels loaded from: {meta_kernel_path}")

# Import classes and methods
from .tides import Body, BodyResponse
from .planet import Earth

# Welcome message
def _welcome_message():
    """Display welcome message on import."""
    print(f"Welcome to tSPICE v{__version__}")
_welcome_message()
del _welcome_message

# Classes and methods to export
__all__ = ['initialize', 'Body', 'BodyResponse', 'Earth', 'get_data']