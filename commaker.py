# -*- coding: utf-8 -*-
"""
Created on Mon Jun  6 09:32:44 2022

@author: ak555
"""
import configparser

import pandas as pd

"""
Functions
---------

template_reader:
    Reads the xyz file and outputs a list of atoms that need to be described
    by the basis sets and the main body of the xyz file.
    
bs_applicator:
    Takes a list of atoms and then specifies basis set functions for each of
    them. It also specifies pseudopotentials for large enough atoms.
    
header_creator:
    This creates the upper portion of the com file, from the supercomputer 
    instructions to the charge and multiplicity.
    
assembler:
    This takes all of the functions from above and concatenates them. 

Notes
-----

These are seperated out into seperate functions for ease of testing and
modification, but this could easily be done as a big function. Maybe one large
function would be better stylistically. 

1. Make template_reader more robust to accept blank lines in the front
2. Allow functionality for a split basis set. 
"""

def template_reader(configFile='config.ini'):
    """

    Parameters
    ----------
    template : STRING, optional
        This will read the template xyz file by default or any xyz file.

    Returns
    -------
    contents : STRING
        This returns the body of the xyz file as a string
    thinned_atoms : LIST
        This returns a list of the atoms included in the xyz file

    """

    config = configparser.ConfigParser()
    config.read(configFile)
    template = config['Files']['template']
    
    read_template = open(template, "rt")
    contents = read_template.read() 
    stripped_contents = contents # I need to have my cake and eat it too.
    read_template.close()
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
    [thinned_atoms.append(x) for x in atomic_species if x not in thinned_atoms]
    ##########

    return contents, thinned_atoms


def bs_applicator(atom_list, configFile='config.ini', pseudo_cutoff=18):
    """

    Parameters
    ----------
    atom_list : LIST
        List of atoms to be included with the basis sets
    basis_set : STRING, optional
        This is the path to the wanted basis set i.e. 'bs/def2-tzvpd.1.gbs' 

    Returns
    -------
    basis_args : STRING
        A string containing the needed basis set arguments with the pseudo-
        potentials added to the end. The basis sets are added in no particular
        order. 

    """
    
    config = configparser.ConfigParser()
    config.read(configFile)
    basis_set = config['Files']['basis_set']
    
    ptable = 'Periodic Table of Elements.csv'
    ptable = pd.read_csv(ptable)
    pseudo_cutoff_ptable = ptable[(ptable['AtomicNumber'] > pseudo_cutoff)]
    PCP_series = pd.Series(pseudo_cutoff_ptable['Symbol']).values
    
    
    atoms = atom_list
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
    
    return basis_args


def header_creator(configFile='config.ini'):
    """

    Parameters
    ----------
    configFile : STRING, optional
        This is the filename of the wanted configuration file of the following
        format:
            
            [Properties]
            charge = 0
            #Total charge of your species in the template
            multiplicity = 1
            #Total multiplicity of your species in the template
            processors = 8
            memory = 16gb
            header = ./headers/default
            title = output test batch

    Returns
    -------
    header : STRING
        This is the top portion of the com file up to the charge/multiplicity

    """
    config = configparser.ConfigParser()
    config.read(configFile)
    # I probably should follow a more strict definition of a header
    procs = '%nprocshared='+config['Properties']['processors']+'\n'
    mem = '%mem='+config['Properties']['memory']+'\n'
    
    header_path = config['Properties']['header']
    read_default_header = open(header_path, "rt")
    instructions = read_default_header.read()+'\n\n'
    read_default_header.close()
    
    title = config['Properties']['title']
    
    charge = config['Properties']['charge']
    multiplicity = config['Properties']['multiplicity']
    charge_multiplicity = charge+" "+multiplicity
    
    header = procs+mem+instructions+title+'\n\n'+charge_multiplicity+'\n'
    
    return header
    

def assembler():
    """

    Returns
    -------
    file : STRING
        This assembles all the parts of the com file created with the other
        functions. 

    """
    xyz, thin_atoms = template_reader()
    basissets = bs_applicator(thin_atoms)
    header = header_creator()
    
    file = header+xyz+basissets+'\n\n' # Gaussian needs blank lines
    
    return file




def main(outputfilename='o.com'):
    text_file = open(outputfilename, "w")
    text_file.write(assembler())
    text_file.close()

if __name__ == '__main__':
    main()