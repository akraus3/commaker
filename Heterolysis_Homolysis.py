import configparser

import pandas as pd
import numpy as np

from Commaker import commaker

"""
Notes
-----

I envision this being run after the main optimization has concluded. 
The first of which is that I am not positive whether or not I should
optimize this or not after the split. In terms of relative energy,
would not the relative energies stay the same? I could then run 
a mere ZPE calculation on the cleaved bits. 

BBB uses strict indexes to identify the atoms that need to be removed.
While that does make it easier programmatically, I would prefer to make
this easier for myself and any other users down the line by using math
and logic to find the atoms that need to be removed. This leads to some
assumptions:

    1. The leaving H is always the furthest from the methane's core

I also don't know if I should optimize the split parts or not. 
-----

Requirements
------------
    
    1. An input .com file
    2. A config.ini file
    3. ./headers/ folder with a text file with Gaussian instructions. IE :
        # opt=(maxcycles=500,noeigentest) freq m06/gen pseudo=read
    4. ./bs/ folder with basis sets downloaded from the exchange.
    5. Commaker.py in the same directory. 
    6. A 'Periodic Table of Elements'.csv file in the root directory. I got
        this from GoodmanSciences at Github.

https://gist.github.com/GoodmanSciences/c2dd862cd38f21b0ad36b8f96b4bf1ee
------------ 
"""

filename = "o.com" # Keep this around just for independent usage.

def Split_Atoms(input_file):
    """

    Parameters
    ----------
    input_file : STRING
        This is the base file. The base file should ideally be produced by 
        commaker, but can be any .com file. 
    
    Returns
    -------
    methyl_fname : STRING
        Returns the filename of the .xyz file produced for the methyl group.
    metalComplex_fname : STRING
        Returns the filename of the .xyz file produced for the metal complex. 
    charge : STRING
        Returns the base file's reported charge. Must be a string of an 
        integer. 
    multiplicity : STRING
        Returns the base file's reported multiplicity. Must be a string of an
        integer. 

    """

    #region load the PTable and base file
    ptable = 'Periodic Table of Elements.csv'
    ptable = pd.read_csv(ptable)
    p_symbols = pd.Series(ptable['Symbol']).values
    metals = ptable[(ptable['Metal'] == 'yes')]
    metals = metals['Symbol']
    read_input = open(input_file, "rt")
    contents = read_input.read()
    read_input.close()
    #endregion

    #region gather xyz_data, charge, and multiplicity
    mol_lines = contents.split('\n') # Splits the string into lines
    xyz_data = []

    # This iterates through the input file to gather the needed information.
    for iii in range(len(mol_lines)):
        # This gathers the charge/multiplicity from the file
        cm_stripper = mol_lines[iii].split(' ')
        if len(cm_stripper) == 2:
            # Don't remember the reasoning, but it doesn't work without this.
            if (len(cm_stripper[0]) + len(cm_stripper[1])) < 3:
                charge = cm_stripper[0]
                multiplicity = cm_stripper[1]
        for symbol in p_symbols: # Get the lines with the atomic symbols
            spaced_symbol = symbol+' '
            if spaced_symbol in mol_lines[iii] and not len(mol_lines[iii]) < 20:
                xyz_coords = mol_lines[iii]
                xyz_data.append(xyz_coords)
    #endregion

    #region format atomic coordinates, find metal core
    metal_coordinates = '' # The coordinates of the metal core
    coordinates = []

    for iii in range(len(xyz_data)):
        atom_coordinates = []
        line = xyz_data[iii].split(' ')
        for iii in range(len(line)):
            if len(line[iii]) > 0: # Remove the extras
                if len(line[iii]) > 2: # If it is one of the coordinates
                    atom_coordinates.append(float(line[iii]))
                else: # if it is just the atomic symbol
                    atom_coordinates.append(line[iii]) 
        
        for metal in metals: # Check for the metal core assumes just one metal
            if atom_coordinates[0] == metal:
                metal_coordinates = atom_coordinates

        coordinates.append(atom_coordinates)
    
    
    del(metal_coordinates[0]) # Remove the atomic symbol for ease of math
    metal_coordinates = np.array(metal_coordinates) # math reasons
    #endregion

    #region find the methyl carbon
    shortest_distance = 2.3 # arbitrary distance to reduce calculations
    leavingCarbon = []
    
    for coordinate in coordinates:
        symbol = coordinate[0]
        xyz_coords = np.array(coordinate[1:])
        if symbol == 'C': # Try to find a Carbon close to the metal core
            x = 0
            difference = metal_coordinates - xyz_coords
            
            distance = np.linalg.norm(difference)
            
            if distance < shortest_distance: # Iterates through each, records
                for H_coord in coordinates: # Next find out if 2 H are close
                    H_coord_symbol = H_coord[0]
                    H_xyz_coords = np.array(H_coord[1:])

                    if H_coord_symbol == 'H':
                        C_H_vector = xyz_coords - H_xyz_coords
                        distance = np.linalg.norm(C_H_vector)

                        if distance < 1.3: #The max distance shown was 1.27
                            x +=1

            if x > 3: #Be absolutely positive that I have the methane
                leavingCarbon = coordinate 
    
    C_coords = np.array(leavingCarbon[1:])
    #endregion

    #region Find the leaving Hydrogen in the transition state. 
    temp_coordinates = []
    temp_distances = []
    for coordinate in coordinates:
        symbol = coordinate[0]
        xyz_coords = np.array(coordinate[1:])

        if symbol == 'H':
            
            difference = xyz_coords - C_coords
            distance = np.linalg.norm(difference)

            if distance < 1.30: # Max bond distance between C and H angstroms
                temp_coordinates.append(coordinate)
                temp_distances.append(distance)

    # Find the coordinates of the leaving Hydrogen
    transition_index = [0, 0]
    for iii in range(len(temp_distances)):
        if temp_distances[iii] > transition_index[0]:
            transition_index = [temp_distances[iii], iii]
    # temp_coordinates will be included in the carbon file, this line just
    # removes the hydrogen that will belong to the metal complex.
    del(temp_coordinates[transition_index[1]])
    
    leavingCarbon = [leavingCarbon]
    for H in temp_coordinates:
        leavingCarbon.append(H)
    #endregion

    #region produce the strings of .xyz files that are needed. 
    methylXYZ = ''
    cmplxXYZ = ''

    for coordinate in xyz_data:
        methyl_coordinate = False # More sorting for the methyl group
        for atom in leavingCarbon: # start filtering through the data again.
            symbol = atom[0]
            x = str(atom[1])
            y = str(atom[2])
            z = str(atom[3])
            
            if (symbol in coordinate and x in coordinate and y in coordinate
                and z in coordinate): # Sort for the methyl group
                methylXYZ += coordinate+'\n'
                methyl_coordinate = True

        if not methyl_coordinate:
            cmplxXYZ += coordinate+'\n'

    # Just some neccessary formatting to work with commaker. 
    methylXYZ += '\n'
    cmplxXYZ += '\n'

    fname = input_file.split('.')[0] # Get rid of the fiel extension
    methyl_fname = fname+"MethylLG.xyz" 
    metalComplex_fname = fname+"CoreComplex.xyz"
    #endregion

    #region write the text files
    text_file = open(methyl_fname, "w")
    text_file.write(methylXYZ)
    text_file.close()

    text_file = open(metalComplex_fname, "w")
    text_file.write(cmplxXYZ)
    text_file.close()
    #endregion

    return methyl_fname, metalComplex_fname, charge, multiplicity

            
def charge_multiplicity(base_charge, base_multiplicity, split_type='homolysis'):
    """

    Parameters
    ----------
    base_charge : STRING
        This is the charge of the input file before the split. Must be an INT in
        a STRING.  
    base_multiplicity : STRING
        This is the multiplicity of the input file before the split. Must be an 
        INT in a STRING. 
    split_type : STRING
        Must be 'homolysis', 'heterolysisNegativeMetal', or 
        'heterolysisPositiveMetal'. 
        TODO ak555 add functionality for the heterolysis.

    Returns
    -------
    metal : TUPLE
        Returns a TUPLE of two STRINGS: the metal's charge and multiplicity in 
        that exact order. 
    carbon : TUPLE
        Returns a TUPLE of two STRINGS: the carbon's charge and multiplicity in 
        that exact order. 

    """
    # Basically stole this from BradenBillyBondFuncs. Makes sense, but I need to
    # finish it. 
    if "homolysis" in split_type:
        metal_multiplicity = int(base_multiplicity) + 1
        carbon_multiplicity = 2 # Should be a doublet. 
        metal_charge = int(base_charge)
        carbon_charge = 0
    else:
        pass


    # TODO ak555 Finish this up. I need to mess with the heterolysis part

    # Messy, and this just digs my pit even deeper. 
    metal = (str(metal_charge), str(metal_multiplicity))
    carbon = (str(carbon_charge), str(carbon_multiplicity))

    return metal, carbon


def main(input, configFile='config.ini'): # Rename this?
    """

    Parameters
    ----------
    input : STRING
        This is the input file. Ideally a .com file produced by commaker, but
        could be any .com file.   
    configFile : STRING, optional
        This is the filename for the config file. The default is 'config.ini'.

    Returns
    -------
    None

    """

    methyl_fname, cmplx_fname, charge, multiplicity = Split_Atoms(input)
    

    methyl_outputfname = methyl_fname.split('.')[0]+'.com'
    cmplx_outputfname = cmplx_fname.split('.')[0]+'.com'
    
    # Read the config file for some of the necessary parameters. 
    config = configparser.ConfigParser()
    config.read(configFile)
    split_type = config['Split']['split_type']

    metal, carbon = charge_multiplicity(charge, multiplicity, split_type)

    basis_set = config['Files']['basis_set']
    title = config['Properties']['title']
    carbon_charge = carbon[0]
    carbon_multiplicity = carbon[1]
    metal_charge = metal[0]
    metal_multiplicity = metal[1]
    procs = '%nprocshared='+config['Properties']['processors']+'\n'
    mem = '%mem='+config['Properties']['memory']+'\n' 
    M_header_path = config['Split']['M_split_header']
    C_header_path = config['Split']['C_split_header']

    # Process each file for commaker. 
    commaker(methyl_fname, basis_set, title, carbon_charge, carbon_multiplicity, 
             procs, mem, C_header_path, methyl_outputfname)
    
    commaker(cmplx_fname, basis_set, title, metal_charge, metal_multiplicity, 
             procs, mem, M_header_path, cmplx_outputfname)
    
    return methyl_outputfname, cmplx_outputfname




if __name__ == '__main__':

    main(filename)