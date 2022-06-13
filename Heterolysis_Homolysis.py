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

    1. The C atom of the methane is always the closest C atom to the metal
    2. The leaving H is always the furthest from the methane's core

These have some clear problems. The first being the most problematic of them
all. 

"""

filename = "o.com"

def Split_Atoms(input_file):

    ptable = 'Periodic Table of Elements.csv'
    ptable = pd.read_csv(ptable)
    p_symbols = pd.Series(ptable['Symbol']).values
    metals = ptable[(ptable['Metal'] == 'yes')]
    metals = metals['Symbol']

    read_input = open(input_file, "rt")
    contents = read_input.read()
    read_input.close()


    mol_lines = contents.split('\n') # Splits the string into lines
    xyz_data = []
    atoms = []
    contents = ''


    for iii in range(len(mol_lines)):
        #This gathers the charge/multiplicity from the file
        cm_stripper = mol_lines[iii].split('  ')
        if len(cm_stripper) == 2:
            charge = cm_stripper[0]
            multiplicity = cm_stripper[1]

        for symbol in p_symbols:
            spaced_symbol = symbol+' '
            if spaced_symbol in mol_lines[iii] and not len(mol_lines[iii]) < 20:
                xyz_coords = mol_lines[iii]
                
                atoms.append(symbol)
                xyz_data.append(xyz_coords)


    metal_coordinates = ''
    coordinates = []


    for iii in range(len(xyz_data)):

        atom_coordinates = []

        line = xyz_data[iii].split(' ')

        for iii in range(len(line)):
            if len(line[iii]) > 0:
                if len(line[iii]) > 2:
                    atom_coordinates.append(float(line[iii]))

                else:
                    atom_coordinates.append(line[iii])


        for metal in metals:
            if atom_coordinates[0] == metal:
                metal_coordinates = atom_coordinates

        coordinates.append(atom_coordinates)
    
    
    del(metal_coordinates[0])
    metal_coordinates = np.array(metal_coordinates)

    shortest_distance = 10
    sd_full = []


    for coordinate in coordinates:
        symbol = coordinate[0]
        xyz_coords = np.array(coordinate[1:])

        if symbol == 'C':
            
            difference = metal_coordinates - xyz_coords

            distance = np.linalg.norm(difference)

            if distance < shortest_distance:
                shortest_distance = distance
                sd_full = coordinate
        

    C_coords = np.array(sd_full[1:])

    temp_coordinates = []
    temp_distances = []

    for coordinate in coordinates:
        symbol = coordinate[0]
        xyz_coords = np.array(coordinate[1:])


        if symbol == 'H':
            
            difference = xyz_coords - C_coords
            distance = np.linalg.norm(difference)

            if distance < 1.30:
                temp_coordinates.append(coordinate)
                temp_distances.append(distance)

    print(temp_coordinates)
    print(temp_distances)

    transition_index = [0, 0]
    for iii in range(len(temp_distances)):
        if temp_distances[iii] > transition_index[0]:
            transition_index = [temp_distances[iii], iii]
    
    del(temp_coordinates[transition_index[1]])
    sd_full.append(temp_coordinates)

    methyl_coords = temp_coordinates

    print(methyl_coords) 



            






Split_Atoms(filename)
