from collections import defaultdict, Counter
import re

def extract_year(id_string):
    """
    Extracts the year from the sequence ID.
    Assumes that the year is a 4-digit number following a pipe '|' in the ID.
    """
    match = re.search(r'[|_](\d{4})\b', id_string)
    if match:
        return int(match.group(1))
    return None

# Function to decide if the group should be kept or removed
def filter_group_by_year(group, min_year=2020, threshold=100):
    # If the group is smaller than 1% of total rows
    if len(group) < threshold:
        years = group['year']
        # Check if years end in min_year or cross min_year
        if years.max() == min_year:
            if years.min() <min_year:
                return False
        if years.max() < min_year:
            return False         
    return True

# def simulate_year_bars(years, total_width=125, start_year=1900, end_year=2025):
#     """
#     Creates a string of dashes representing the year range based on a list of years.

#     Parameters:
#     - years: List of years.
#     - total_width: Number of characters to represent the year range (default: 40).
#     - start_year: Start of the year range (default: 1985).
#     - end_year: End of the year range (default: 2025).

#     Returns:
#     - bar: A string representing the year span using dashes.
#     """
#     if not years:
#         return ' ' * total_width

#     min_year, max_year = min(years), max(years)
#     start_pos = int(((min_year - start_year) / (end_year - start_year)) * total_width)
#     end_pos = int(((max_year - start_year) / (end_year - start_year)) * total_width)

#     bar = [' ' for _ in range(total_width)]
#     for i in range(start_pos, end_pos + 1):
#         bar[i] = '-'
    
#     return ''.join(bar)

def simulate_year_bars(years, total_width=125, start_year=1900, end_year=2025):
    """
    Creates a string of dashes representing individual years within a given range.

    Parameters:
    - years: List of years.
    - total_width: Number of characters to represent the year range (default: 125).
    - start_year: Start of the year range (default: 1900).
    - end_year: End of the year range (default: 2025).

    Returns:
    - bar: A string representing the years using dashes.
    """
    if not years:
        return ' ' * total_width

    bar = [' ' for _ in range(total_width)]
    for year in years:
        if start_year <= year <= end_year:
            pos = int(((year - start_year) / (end_year - start_year)) * total_width)
            bar[pos] = '-'
    
    return ''.join(bar)


def get_subsequences_with_years(df_samples, kmer_column):
    """
    Extracts subsequences and their associated years from the given DataFrame.

    Parameters:
    - df_samples: pandas DataFrame with columns ['ID', kmer_column].
    - kmer_column: Name of the column containing kmers.

    Returns:
    - subseq_to_years: Dictionary mapping subsequences to associated years.
    - subseq_freq: Counter object for the frequency of subsequences.
    """
    df_clean = df_samples[['ID', kmer_column]].dropna()

    subseq_to_years = defaultdict(list)
    for idx, row in df_clean.iterrows():
        year = extract_year(row['ID'])
        if year:
            subseq_to_years[row[kmer_column]].append(year)
    
    subseq_freq = Counter(df_clean[kmer_column])
    return subseq_to_years, subseq_freq

def get_sorted_subsequences(subseq_freq):
    """
    Sorts subsequences by frequency in descending order.

    Parameters:
    - subseq_freq: Counter object with subsequence frequencies.

    Returns:
    - sorted_subsequences: List of sorted subsequences with their frequencies.
    """
    return sorted(subseq_freq.items(), key=lambda x: x[1], reverse=True)

