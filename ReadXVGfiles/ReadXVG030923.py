# Define a function to read xvg files
import pandas as pd

def read_xvg(file_path):
    data = []

    with open(file_path, 'r') as f:
        for line in f:
            if not line.startswith("#") and not line.startswith('@'):
                new_list = [elem for elem in line.split()]
                data.append(new_list)
    # Convert the list of lists into a Pandas DataFrame
    df = pd.DataFrame(data)
    return df

