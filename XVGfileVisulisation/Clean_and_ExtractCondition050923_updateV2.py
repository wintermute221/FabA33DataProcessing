"""
Script Name: Index writing Script
Author: Yuhan Wang
Creation Date: 4th Sep, 2023

Description:
This script aims to extract the data info from a large number of XVG files and then merge them together wihtin one dataframe with the condition as the column title

- This script is intended for educational purposes only.
- Use at your own risk.

"""

# clean & extract conditions
##Steps:
##1) get rid of the suffix .xvg 2) get rid of the replica naming r2 3) get rid of rmsd 4) Get rid of ##the loop if there is any 5) get rid of any left letters like M or K


# define a function to clean the filename so as to get the 54 condition names

import re

pattern_loop = r'loop\d{1}'  # Pattern to remove 'loop' followed by a digit
pattern_rmsd = r'rmsd'  # Pattern to remove 'rmsd'


def clean_filename(filename):
    # Remove the '.xvg' extension
    filename = filename.replace('.xvg', '')

    # get rid of i.g. _r2
    filename = re.sub(r'r\d{1}', '', filename)

    # Remove 'loop' followed by a digit
    filename = re.sub(pattern_loop, '', filename, flags=re.IGNORECASE)

    # Remove 'rmsd'
    filename = re.sub(pattern_rmsd, '', filename, flags=re.IGNORECASE)

    # Remove Replicas
    filename = re.sub(r'Replicas', '', filename)

    # Remove protonation_pH3_5 & protonation_pH4_5
    filename = re.sub(r'protonation_pH3_5', '', filename)
    filename = re.sub(r'protonation_pH4_5', '', filename)

    # Remove any non-numeric characters except underscores
    # filename = re.sub(r"[^0-9_]", "", filename)

    # Remove leading and trailing underscores
    filename = filename.strip('_')

    return filename