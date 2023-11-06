"""
Script Name: Native contact Script
Author: Yuhan Wang
Creation Date: 6th Nov, 2023
Updated version

Description:
This script aims to calculate the native contacts across the trajectory

Calculate contacts for all the alpha-carbons in the protein
Define the contact radius cutoff at 8 Angstrom
The contacts.q1q2 function uses the radius_cut_q method to calculate the fraction of native contacts for a conformation by determining that atoms i and j are in contact if they are within a given radius

Reference: https://userguide.mdanalysis.org/1.0.1/examples/analysis/distances_and_contacts/contacts_q1q2.html

- This script is intended for educational purposes only.
- Use at your own risk.

"""

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import pandas as pd
import os

# Change to the targeting directory

os.chdir(
    "/Users/wintermute/Library/CloudStorage/GoogleDrive-wintermute.backup@gmail.com/My Drive/MacbookPro/GoogleDrive/UCL/PaulGroup/MolecularDynamics/Gromacs/FabA33/DataAnalysisAll200323/NativeContacts250423/TryDifferentCuttingMethods051123/UnproblematicMatching061123/54DCD140623")

# define directory
path = os.getcwd()

# read gro file
pdb = "Fab.pdb"

# create a list for saving average native contacts
data = []
cols = ['Condition', 'Average native contact (last 50ns)']

# loop over every dcd file in the directory
for file in os.listdir(path):

    # select the file that ends with .dcd
    # reference: https://stackoverflow.com/questions/10377998/how-can-i-iterate-over-files-in-a-given-directory
    if file.endswith(".dcd"):

        try:

            # read dcd trajecotry
            dcd = file
            u = mda.Universe(pdb, dcd)

            # get the filename without the extension
            filename = os.path.splitext(file)[0]

            q1q2 = contacts.q1q2(u, 'name CA', radius=8).run()

            q1q2_df = pd.DataFrame(q1q2.timeseries,
                                   columns=['Frame',
                                            'Q1',
                                            'Q2'])

            # plot native contacts
            q1q2_df.plot(x='Frame')
            plt.ylabel('Fraction of native contacts')

            plt.savefig(filename + 'q1q2.png')

            # Calculate the average native contacts of the last 500 frames
            average_contacts = q1q2_df.tail(500).mean()

            # Append the Q1 values into the dataframe
            data.append([filename, average_contacts["Q1"]])


        except:

            pass

        continue

    else:
        continue

    # Save the average last 50ns native contacts into an excel file

    # df = pd.DataFrame(data, columns=cols)
    # df.to_excel('NativeContacts_AlphaCarbon061123.xlsx')