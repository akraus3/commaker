import os
import subprocess
import configparser

import Heterolysis_Homolysis
from Commaker import commaker


def main(configFile='config.ini'):

    #region commaker input arguments
    config = configparser.ConfigParser()
    config.read(configFile)
    basis_set = config['Files']['basis_set']
    title = config['Properties']['title']
    charge = config['Properties']['charge']
    multiplicity = config['Properties']['multiplicity']    
    procs = '%nprocshared='+config['Properties']['processors']+'\n'
    mem = '%mem='+config['Properties']['memory']+'\n' 
    header_path = config['Properties']['header']
    #endregion

    # Much of this is stolen from BatchBadenBond because I liked its output.
    if not os.path.isdir('results'):
        os.mkdir('results')

    for file in os.listdir('inputs/'):
        
        full_filename = 'inputs/'+file
        filename = file.split('.')[0]

        if not os.path.isdir('results/'+filename):
            os.mkdir('results/'+filename)
        com_filename = 'results/'+filename+'/'+filename+'.com'
        
        commaker(full_filename, basis_set, title, charge, multiplicity, procs, mem,
                 header_path, com_filename)
        
        C_fname, metal_fname = Heterolysis_Homolysis.main(com_filename)

        subprocess.check_call('rung16 %s' % (com_filename), shell=True)
        subprocess.check_call('rung16 %s' % (C_fname), shell=True)
        subprocess.check_call('rung16 %s' % (metal_fname), shell=True)
    
        

        
        

        

main()