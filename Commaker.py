# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 09:32:44 2022

@author: ak555
"""
import configparser

import pandas as pd

"""
Notes
-----

Generally, I see this program as being used in conjunction with a config file
of the following format:
    
    [Properties]
    charge = 0
    #Total charge of your species in the template
    multiplicity = 1
    #Total multiplicity of your species in the template
    processors = 8
    memory = 16gb
    header = ./headers/default
    title = output test batch
    [Files]
    template = template.mol
    basis_set = bs/def2-tzvpd.1.gbs

However, it has been seperated for improved usability, because if I want to 
automate the calculations of hundreds of atoms, I need to be able to run this
with many different configurations for output filenames and such. 

I just added .gjf file functionality to this based on the files openbabel 
would output given the MASON xtb optimized structures and the command:
    obabel *.mol -ogjf -m
This is great because it solved my charge and multiplicity problem, or, at
the very least, it pushes it down the road.
-----

Requirements
------------
    
    1. A template file of an .xyz or .mol or .gjf type
    2. A config.ini file
    3. ./headers/ folder with a text file with Gaussian instructions. IE :
        # opt=(maxcycles=500,noeigentest) freq m06/gen pseudo=read
    4. ./bs/ folder with basis sets downloaded from the exchange.
    5. A 'Periodic Table of Elements'.csv file in the root directory. I got
        this from GoodmanSciences at Github.

https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee
------------ 

MISC
----
TODO AK555
1. Make template_reading func more robust to accept blank lines in the front
2. Allow functionality for a split basis set. 
3. Custom SMD solvent ability?
"""

def commaker(template, basis_set, title, charge, multiplicity, procs, 
             mem, header_path, outputfilename, pseudo_cutoff=18):
    """

    Parameters
    ----------
    template : STRING
        This is the file that will be transformed. 
    basis_set : STRING
        This is the desired basis set for the calculation. 
    title : STRING
        This is the title line required by Gaussian. 
    charge : STRING
        See multiplicity
    multiplicity : STRING
        TODO AK555 : for the atom splitting the charge/multiplicity will change
        Maybe I should make this an int?
    procs : STRING
        The amount of required processors. The format is very particular:
            
            "%mem=16gb\n"
            
    mem : STRING
        The amount of memory needed for the job. The format is very particular:
            
            "%nprocshared=8\n"
            
    header_path : STRING
        I made this so that I could have several different gaussian instruction
        files made that could be used. IE a Transition State header, or one
        with a different temperature. 
    outputfilename : STRING
        This is the path and filename for the .com file that is produced. 
    pseudo_cutoff : STRING, optional
        This is the smallest atomic number that will be given pseudopotentials
        by the program. The default is 18.

    Returns
    -------
    None.

    """

    # GoodmanSciences csv
    ptable = 'Periodic Table of Elements.csv'
    ptable = pd.read_csv(ptable)
    p_symbols = pd.Series(ptable['Symbol']).values
    pseudo_cutoff_ptable = ptable[(ptable['AtomicNumber'] > pseudo_cutoff)]
    PCP_series = pd.Series(pseudo_cutoff_ptable['Symbol']).values   
    
    CM_passed = False

    if not '_io.StringIO' in str(type(template)):
        read_template = open(template, "rt")
        contents = read_template.read() 
        stripped_contents = contents # I need to have my cake and eat it too.
        read_template.close()
    
    # This is the pathway for xyz files of a very clean type
    if template.split('.')[-1] == 'xyz' or '_io.StringIO' in str(type(template)):
        ##########
        # These lines take the template and boil down the template xyz file
        # into just a list of atoms without duplicates
        ##########
        # Remove the coordinates and spaces
        bad_char = ['1', '2', '3', '4', '5', '6', '7',
                    '8', '9', '0', '.', ' ', '-', '\t']
        for char in bad_char:
            stripped_contents = stripped_contents.replace(char, '')
        atomic_species = stripped_contents.split("\n")
        ##########
        # Remove blank lines
        iii = 0
        while iii != len(atomic_species):
            if len(atomic_species[iii]) == 0:
                del(atomic_species[iii])
                iii -= 1
            iii += 1
        ##########
        # Remove duplicate atoms using a list comprehension
        thinned_atoms = []
        atoms = atomic_species # Just a lame fix to keep the line short enough
        [thinned_atoms.append(x) for x in atoms if x not in thinned_atoms]
        ##########
    
    # This is the pathway for the mol files, particularly those created by 
    # MASON.
    elif template.split('.')[-1] == 'mol':
        mol_lines = contents.split('\n') # Splits the string into lines
        xyz_data = []
        atoms = []
        contents = ''
        
        ##########
        # This iterates through each line of the file and sees if it has a real 
        # atomic number that is spaced correctly.
        ##########
        for iii in range(len(mol_lines)):
            for symbol in p_symbols:
                spaced_symbol = ' '+symbol+' '
                if spaced_symbol in mol_lines[iii]:
                    
                    coords = mol_lines[iii].split(spaced_symbol)[0]
                    xyz_coords = symbol+coords # Reorders the coordinates
                    
                    atoms.append(symbol)
                    xyz_data.append(xyz_coords)
        
        # Generates an xyz-esque string
        for iii in range(len(xyz_data)):
            contents += xyz_data[iii]+'\n'
        contents += '\n'
        ##########
        # Remove duplicate atoms using a list comprehension
        thinned_atoms = []
        [thinned_atoms.append(x) for x in atoms if x not in thinned_atoms]
        ##########        
    
    elif template.split('.')[-1] == 'gjf':
        mol_lines = contents.split('\n') # Splits the string into lines
        xyz_data = []
        atoms = []
        contents = ''

        ##########
        # This iterates through each line of the file and sees if it has a real 
        # atomic number that is spaced correctly.
        ##########
        for iii in range(len(mol_lines)):
            #This gathers the charge/multiplicity from the file
            cm_stripper = mol_lines[iii].split('  ')
            if len(cm_stripper[0]) == 3 and not CM_passed:
                cm_stripper = cm_stripper[0].split(' ')
                print(cm_stripper)
                charge = cm_stripper[0]
                multiplicity = cm_stripper[1]
                CM_passed = True

            for symbol in p_symbols:
                spaced_symbol = symbol+' '
                if spaced_symbol in mol_lines[iii]:
                    xyz_coords = mol_lines[iii]
                    
                    atoms.append(symbol)
                    xyz_data.append(xyz_coords)

        # Generates an xyz-esque string
        for iii in range(len(xyz_data)):
            contents += xyz_data[iii]+'\n'
        contents += '\n'
        ##########
        # Remove duplicate atoms using a list comprehension
        thinned_atoms = []
        [thinned_atoms.append(x) for x in atoms if x not in thinned_atoms]
        ##########         

    print(charge)
    print(multiplicity)


    ##########
    # The Basis Set section. 
    atoms = thinned_atoms
    # Open the basis set    
    read_basis = open(basis_set, "rt")
    bs = read_basis.read()
    read_basis.close()
    ##########
    # Sort the text file into a list of three elements
    full_file = bs.split('\n\n')
    label = full_file[0]
    basis = full_file[1]
    pseudopotentials = full_file[2].split("    0")
    ##########
    # Split up the basis set section
    basis_functions = basis.split('\n****')
    ##########
    # Gather the necessary Basis functions based on the needed list of atoms
    needed_functions = []
    for atom in atoms:
        atom = atom+"     0"
        for function in basis_functions:
            if atom in function:
                footer = "\n****"
                corrected_function = function+footer
                needed_functions.append(corrected_function)
    ##########
    # Gather the necessary pseudopotentials based on the list of atoms
    needed_pseudopotentials = []
    bigenough_atoms = []
    for atom in atoms:
        if atom in PCP_series:
            bigenough_atoms.append(atom)  
    # This then gathers the needed pseudopotentials        
    if len(bigenough_atoms) > 0:
        for atom in bigenough_atoms:
            atom = atom.upper()+" " 
            for iii in range(len(pseudopotentials)):
                if atom in pseudopotentials[iii]:
                    header = "\n\n"+atom+"    0"
                    pp = pseudopotentials[iii+1]
                    # The way I sorted the pseudopotentials uses "     0"
                    # as a delimiter which makes a list that is weird
                    pp = pp[:pp.rfind('\n')]
                    pp = header+pp
                    needed_pseudopotentials.append(pp)
                    # All the rest is just cleanup
    ##########
    # Finally, just append to one great big string
    basis_args = label
    for func in needed_functions:
        basis_args = basis_args+func
    for potential in needed_pseudopotentials:
        basis_args = basis_args+potential    
    

    # I probably should follow a more strict definition of a header
    
    read_default_header = open(header_path, "rt")
    instructions = read_default_header.read()+'\n\n'
    read_default_header.close()
    

    charge_multiplicity = charge+" "+multiplicity
    
    header = procs+mem+instructions+title+'\n\n'+charge_multiplicity+'\n'
    
    
    xyz = contents
    basissets = basis_args
    
    file = header+xyz+basissets+'\n\n' # Gaussian needs blank lines
    
    text_file = open(outputfilename, "w")
    text_file.write(file)
    text_file.close() 




def main(outputfilename='o.com', configFile='config.ini'):
    
    config = configparser.ConfigParser()
    config.read(configFile)
    template = config['Files']['template']
    basis_set = config['Files']['basis_set']
    title = config['Properties']['title']
    charge = config['Properties']['charge']
    multiplicity = config['Properties']['multiplicity']    
    procs = '%nprocshared='+config['Properties']['processors']+'\n'
    mem = '%mem='+config['Properties']['memory']+'\n' 
    header_path = config['Properties']['header']
    
    commaker(template, basis_set, title, charge, multiplicity, procs, mem,
             header_path, outputfilename)
    


if __name__ == '__main__':
    main()