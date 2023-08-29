"""
Script Name: Native contact Script
Author: Yuhan Wang
Creation Date: 14th June, 2023

Description:
This script aims to calculate the native contacts across the trajectory

- This script is intended for educational purposes only.
- Use at your own risk.

"""

import numpy as np
import matplotlib.pyplot as plt
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import pandas as pd
import os

# Change to the 54 DCD files' directory

os.chdir(
    '/Users/wintermute/Library/CloudStorage/GoogleDrive-wintermute.backup@gmail.com/My Drive/MacbookPro/GoogleDrive/UCL/PaulGroup/MolecularDynamics/Gromacs/FabA33/DataAnalysisAll200323/NativeContacts250423/NativeContact54runs140623/54DCD140623')

# define directory
path = os.getcwd()

# read pdb
pdb = "Fab.pdb"

# create a list for saving average native contacts
data = []
cols = ['Condition', 'Average native contact (last 50ns)']

# loop over every dcd file in the directory
for file in os.listdir(path):

    # select the file that ends with .dcd
    # reference: https://stackoverflow.com/questions/10377998/how-can-i-iterate-over-files-in-a-given-directory
    if file.endswith(".dcd"):

        # read dcd trajecotry
        dcd = file
        u = mda.Universe(pdb, dcd)

        # get the filename without the extension
        filename = os.path.splitext(file)[0]

        # crude definition of salt bridges as contacts between NH/NZ in ARG/LYS and
        # OE*/OD* in ASP/GLU. You might want to think a little bit harder about the
        # problem before using this for real work.

        # reference: https://docs.mdanalysis.org/1.0.1/documentation_pages/analysis/contacts.html

        # As the comments above, perhaps need to double check whether use this version of defination of salt bridges

        # this defination doesn't work in this case

        sel_basic = "(resname ARG or resname LYS) and (name NH* NZ)"
        sel_acidic = "(resname ASP or resname GLU) and (name OE* OD*)"

        # the (acidic, basic) selections from u, which are assigned from the first frame

        basic = u.select_atoms(sel_basic)
        acidic = u.select_atoms(sel_acidic)

        # radius set as 5.5 A
        ca2 = contacts.Contacts(u, select=(sel_acidic, sel_basic),
                                refgroup=(acidic, basic),
                                radius=5.5, # change raidus here
                                method='radius_cut').run()

        ca2_df = pd.DataFrame(ca2.timeseries,
                              columns=['Frame', 'Contacts from first frame'])

        # print number of average contacts for each condition
        average_contacts = np.mean(ca2.timeseries[501:, 1])

        # append the values into a list first and then convert it to a dataframe
        # this will save a lot of performance time
        # reference: https://stackoverflow.com/questions/31674557/how-to-append-rows-in-a-pandas-dataframe-in-a-for-loop
        data.append([filename, average_contacts])

        # plot native contacts
        ca2_df.plot(x='Frame')
        plt.ylabel('Fraction of contacts')
        plt.title(filename)

        plt.savefig(filename + '.png')

        continue

    else:
        continue

#save the average native contacts for each condition into an excel file
df1 = pd.DataFrame(data, columns=cols)
df1.to_excel('54DCDTest140623.xlsx')