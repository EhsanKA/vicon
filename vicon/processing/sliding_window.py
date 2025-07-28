import numpy as np


def sliding_window_matches(seq1, seq2, window_size=150, threshold=145):
    """
    Generates a binary list indicating whether the number of matches in a sliding window
    exceeds a given threshold.

    Parameters:
        seq1 (str): Reference sequence.
        seq2 (str): Sample sequence.
        window_size (int): Size of the sliding window.
        threshold (int): Threshold for match counts.

    Returns:
        list: A binary list (0 or 1) indicating if the match count exceeds the threshold.
    """
    def encode(sequence):
        mapping = {"A": 0, "T": 1, "C": 2, "G": 3}
        return np.array([mapping.get(char, -1) for char in sequence], dtype=int)  # Map invalid chars to -1

    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length.")

    encoded_seq1 = encode(seq1)
    encoded_seq2 = encode(seq2)

    # Initialize match count for the first window
    current_match_count = np.sum(
        (encoded_seq1[:window_size] == encoded_seq2[:window_size]) & (encoded_seq1[:window_size] != -1)
    )
    match_counts = [current_match_count]

    # Slide the window
    for i in range(window_size, len(seq1)):
        # Subtract the outgoing element
        if encoded_seq1[i - window_size] == encoded_seq2[i - window_size] and encoded_seq1[i - window_size] != -1:
            current_match_count -= 1
        # Add the incoming element
        if encoded_seq1[i] == encoded_seq2[i] and encoded_seq1[i] != -1:
            current_match_count += 1
        match_counts.append(current_match_count)

    # Convert match counts to binary results based on the threshold
    binary_result = (np.array(match_counts) >= threshold).astype(np.uint8)
    return binary_result



def calculate_validity_mask(seq1, seq2, window_size=150):
    """
    Calculates a validity mask for sliding windows, indicating whether each window
    contains only valid characters.

    Parameters:
        seq1 (str): Reference sequence.
        seq2 (str): Sample sequence.
        window_size (int): Size of the sliding window.

    Returns:
        list: A binary mask (0 or 1) where 1 indicates the window is valid (no invalid characters),
              and 0 indicates the window is invalid.
    """
    def encode(sequence):
        mapping = {"A": 0, "T": 1, "C": 2, "G": 3}
        return np.array([mapping.get(char, -1) for char in sequence], dtype=int)  # Map invalid chars to -1

    if len(seq1) != len(seq2):
        raise ValueError("Sequences must have the same length.")

    encoded_seq1 = encode(seq1)
    encoded_seq2 = encode(seq2)

    # Initialize variables
    current_valid_count = 0

    # Calculate initial validity for the first window
    for i in range(window_size):
        if encoded_seq1[i] != -1 and encoded_seq2[i] != -1:  # Both characters are valid
            current_valid_count += 1

    # Initialize mask
    validity_mask = [1 if current_valid_count == window_size else 0]

    # Slide the window
    for i in range(window_size, len(seq1)):
        # Update the outgoing character
        if encoded_seq1[i - window_size] != -1 and encoded_seq2[i - window_size] != -1:
            current_valid_count -= 1

        # Update the incoming character
        if encoded_seq1[i] != -1 and encoded_seq2[i] != -1:
            current_valid_count += 1

        # Append validity result for the current window
        validity_mask.append(1 if current_valid_count == window_size else 0)

    return validity_mask
