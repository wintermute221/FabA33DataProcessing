"""
Script Name: Salt bridge occurrence through trajectories calculation
Author: Yuhan Wang
Creation Date: 27th Nov, 2023

Description:
Workflow:
1. Define a function to calculate salt bridge percentage occurrence based on the dcd and gro input files
2. Use the definition to calculate the percentage occurrence for all the 324 files and save all the data into a dictionary whose key is the dcd file name (condition name)
3. Loop over the dictionary to average the replicas (to note, this step, technically speaking, can be combined into step 2 but based on my experiment it just took too long the calcuation time and very easily induce errors. That is why I put this step here)

The essential ideas here are i) firstly to find all the salt bridges by definition that as long as the distance between N and O of a pair is smaller than 3.2 in any of the frame across the trajectory then it is a salt bridge, ii) Then we are interested in how which salt bridge last in the final 500 frames and which one broke. The pair that has 0 occurrence percentage in the final excel file (/Users/wintermute/Library/CloudStorage/GoogleDrive-wintermute.backup@gmail.com/My Drive/MacbookPro/GoogleDrive/UCL/PaulGroup/MolecularDynamics/Gromacs/FabA33/DataAnalysisAll200323/SaltBridge261123/SaltBridges54Conditions271123.xlsx) means it was a salt bridge because in that condition at least one frame (out of 1000) has the distance smaller than 3.2, but it has no frame in the last 500 ones that has a distance smaller than 3.2



- This script is intended for educational purposes only.

"""


import os
import re
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from glob import glob
import numpy as np
import matplotlib.pyplot as plt


def calculate_salt_bridge_percentage_occurrence(u, group_a, group_b, radius=3.2, last_n_rows=500):
    """
    Calculate salt bridge distances and percentage occurrence of non-NaN values in the last N rows.

    Parameters:
    - u: MDAnalysis Universe
    - group_a: First group of atoms
    - group_b: Second group of atoms
    - radius: Cutoff radius for salt bridge calculation
    - last_n_rows: Number of last rows to consider for percentage occurrence calculation

    Returns:
    - result_df: DataFrame containing salt bridge distances
    - percentage_occurrence: Series with the percentage occurrence of non-NaN values in each column
    """
    
    # Function to calculate salt bridge distances
    def contacts_within_cutoff(u, group_a, group_b, radius=3.2):
        timeseries = []
        pos = [[], []]
        for ts in u.trajectory:
            # calculate distances between group_a and group_b
            dist = contacts.distance_array(group_a.positions, group_b.positions)
            # determine which distances <= radius and save all the qualified salt bridges across frames
            for i, lin in enumerate(dist):
                for j, col in enumerate(lin):
                    if col < radius and col != 0.0:
                        pos[0].append(i) # acids
                        pos[1].append(j) # basics
                        # group variables into one string
                        result_string = group_a[i].resname + str(group_a[i].resid) + '-' + group_b[j].resname + str(group_b[j].resid)
                        element = [ts.frame, result_string, round(col, 2)]
                        timeseries.append(element)
        return np.array(timeseries)

    # Calculate salt bridge distances and create DataFrame
    ca = contacts_within_cutoff(u, group_a, group_b, radius)
    ca_df = pd.DataFrame(ca)
    ca_df.columns = ['frame', 'salt bridge', 'distance']
    
    # Extract unique values for the new DataFrame
    unique_salt_bridges = ca_df['salt bridge'].unique()
    unique_frames = ca_df['frame'].unique()
    
    # Create a new DataFrame with unique values as columns and index
    new_df = pd.DataFrame(index=unique_frames, columns=unique_salt_bridges)
    
    # Fill in the distances based on matching frame and salt bridge
    for _, row in ca_df.iterrows():
        new_df.at[row['frame'], row['salt bridge']] = row['distance']

    # Select the last 500 rows
    last_n_rows_df = new_df.tail(last_n_rows)
    
    # Calculate the occurrence of non-NaN values in each column
    occurrence_counts = last_n_rows_df.count()

    # Calculate the total number of rows
    total_rows = len(last_n_rows_df)

    # Calculate the percentage of non-NaN values
    percentage_occurrence = (occurrence_counts / total_rows)

    return new_df, percentage_occurrence

# Example usage:
# result_df, percentage_occurrence = calculate_salt_bridge_percentage_occurrence(u, acidic, basic, radius=3.2, last_n_rows=500)
# percentage_occurrence

# Step 1: Calculate percentage occurrence for all files

# change to your own working directory
os.chdir("/Users/wintermute/Library/CloudStorage/GoogleDrive-wintermute.backup@gmail.com/My Drive/MacbookPro/GoogleDrive/UCL/PaulGroup/MolecularDynamics/Gromacs/FabA33/DataAnalysisAll200323/SaltBridge261123/DCD54runs/UnproblematicMatching061123")
path = os.getcwd()

fab = "Fab.pdb"
d_salt_bridge = {}  # dictionary that will hold them 

for file in glob(os.path.join(path, "*.dcd")):
    try:
        dcd = file
        u = mda.Universe(fab, dcd)
        filename = os.path.splitext(os.path.basename(file))[0]

        sel_basic = "(resname ARG or resname LYS) and (name NH* NZ)"
        sel_acidic = "(resname ASP or resname GLU) and (name OE* OD*)"

        basic = u.select_atoms(sel_basic)
        acidic = u.select_atoms(sel_acidic)

        _, percentage_occurrence = calculate_salt_bridge_percentage_occurrence(u, acidic, basic, radius=3.2, last_n_rows=500)
        d_salt_bridge[filename] = percentage_occurrence

    except Exception as e:
        print(f"Error processing {file}: {e}")

# Step 2: Loop over the dictionary to average the replicas
d_salt_bridge_54condition = {}

run54 = pd.read_excel("/Users/wintermute/Library/CloudStorage/GoogleDrive-wintermute.backup@gmail.com/My Drive/MacbookPro/GoogleDrive/UCL/PaulGroup/MolecularDynamics/Gromacs/FabA33/DataAnalysisAll200323/SaltBridge261123/54ConditionName120923.xlsx", index_col=0)
run54 = run54['Condition']

for condition in run54:
    replicas = [d_salt_bridge[key] for key in d_salt_bridge if re.search(condition, key)]
    if replicas:
        d_salt_bridge_54condition[condition] = pd.DataFrame(replicas).mean(axis=0)

# Filter out empty values
d_salt_bridge_44condition = {key: value for key, value in d_salt_bridge_54condition.items() if not value.empty}


# # Step 3: Save to Excel
# with pd.ExcelWriter('SaltBridgesMatching.xlsx') as writer:
#     # Iterate over the dictionary
#     for sheet_name, series in d_salt_bridge_54condition.items():
#         # Convert the series to a DataFrame and reset the index
#         df = pd.DataFrame({'Salt Bridge': series.index, 'Percentage of occurrence in the last 50ns': series.values})
#         df.to_excel(writer, sheet_name=sheet_name, index=False)
