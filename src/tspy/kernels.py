import spiceypy as spy
import json
import os

#Function to download the kernels
def download_kernels(data_dir):

    #Directory for the kernels and the configuration file
    kernel_dir = os.path.join(data_dir, 'kernels')
    kernel_config = os.path.join(data_dir, 'spice_kernels.json')

    #print(f'current_dir: {current_dir}')
    #print(f'data_dir: {data_dir}')
    #print(f'kernel_dir: {kernel_dir}')
    #print(f'kernel_config: {kernel_config}')

    #Read the JSON file
    with open(kernel_config, 'r') as f:
        config = json.load(f) 	#Reads the file and turns it into a Python dictionary!

    #If the kernel directory does not exist, create it and download the kernels
    if not os.path.exists(kernel_dir):
        os.makedirs(kernel_dir)
        print(f"Directory for Kernels created.")

        #Start downloading the kernel files
        for i,k in enumerate(config['kernels']):
            try:
                print(f"Downloading the file {k['filename']} from {k['url']}...\n")
                os.system(f"wget -P {kernel_dir} {k['url']}")
                print(f"The file {k['filename']} has been downloaded.\n")
            except:
                print(f"Error downloading the file {k['filename']} from {k['url']}.\n")

    #If the directory exists, verify that the files are there already
    else:
        print(f"Directory for Kernels already exists.")
        
        #Verify that the files are there already or download them if they are not
        for i,k in enumerate(config['kernels']):
            if k['filename'] in os.listdir(kernel_dir):
                print(f"The file {k['filename']} already exists in the kernel directory.\n")
            else:
                #Download the kernel file if it does not exist
                try:
                    print(f"Downloading the file {k['filename']} from {k['url']}...\n")
                    os.system(f"wget -P {kernel_dir} {k['url']}")
                    print(f"The file {k['filename']} has been downloaded.\n")
                except:
                    print(f"Error downloading the file {k['filename']} from {k['url']}.\n")

#Function to write the meta kernel
def write_meta_kernel(data_dir):

    #Create the path to the meta kernel
    meta_kernel_path = os.path.join(data_dir, 'meta_kernel')

    #Create the path to the kernel directory
    kernel_dir = os.path.join(data_dir, 'kernels')

    #If the meta kernel does not exist, create it
    if not os.path.exists(meta_kernel_path):

        with open(meta_kernel_path, 'w') as f:

            f.write('KPL/MK\n')	#Type of file
            f.write('\n')	#Beginning of data
            f.write('\\begindata\n')	#Beginning of data

            f.write(f"PATH_VALUES = ('{kernel_dir}')\n")	#Path to the kernels
            f.write("PATH_SYMBOLS = ('KERNELS')\n")	#Symbol to indicate the path
            f.write('\n')	#New line	

            #Writing the kernels
            f.write('KERNELS_TO_LOAD = (\n')
            for k in os.listdir(kernel_dir):
                f.write(f"'$KERNELS/{k}',\n")
            f.write(')\n')
        print(f"Meta kernel created at {meta_kernel_path}.")

    #If the meta kernel already exists, return its path
    else:
        print(f"Meta kernel already exists at {meta_kernel_path}.")

    #It always returns the path to the meta kernel
    return meta_kernel_path
