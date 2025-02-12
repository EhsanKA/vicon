import pandas as pd

def read_wig(wig_file_path):
    """
    Reads a WIG file and returns the values as a pandas Series with positions as index.
    """
    positions = []
    values = []
    
    print(f"Reading WIG file: {wig_file_path}")
    with open(wig_file_path, 'r') as wig_file:
        for line in wig_file:
            if line.startswith('variableStep'):
                continue
            parts = line.strip().split()
            if len(parts) == 2:
                position = int(parts[0])
                value = float(parts[1])
                positions.append(position)
                values.append(value)
    data_series = pd.Series(data=values, index=positions)
    print(f"Read {len(data_series)} values from {wig_file_path}")
    return data_series