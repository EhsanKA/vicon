import matplotlib.pyplot as plt

def plot_non_gap_counts(df_counts, title='Non-Gap Character Counts at Each Position', save="", logger=None):
    """
    Plots the counts of non-gap characters at each position and saves the plot based on the title.

    Parameters:
    df_counts (pd.DataFrame): DataFrame containing counts of non-gap characters at each position.
    title (str): Title of the plot.
    """
    plt.figure(figsize=(12, 6))
    plt.bar(df_counts.index, df_counts['NonGapCount'], color='blue')
    plt.xlabel('Position')
    plt.ylabel('Number of Non-Gap Characters')
    plt.title(title)
    plt.tight_layout()
    
    # Generate a filename based on the title
    # filename = title.replace(' ', '_') + '.png'
    
    # Save the plot to a file
    plt.savefig(save)
    if logger:
        logger.info(f"Plot saved as {save}")
    else:
        print(f"Plot saved as {save}")
    
    # Close the plot to free up memory
    plt.close()


import matplotlib.pyplot as plt
import numpy as np

def plot_rel_cons(df, kmer_size=150, threshold=3, save_path='.', sample_name='sample'):
    
    # Assuming df has columns 'position', 'value_relaxed', and 'value_other'
    x = df.columns.values.flatten()
    y = df.sum().values.flatten()
    # Finding change indices for y
    change_indices_y = np.where(y[:-1] != y[1:])[0]
    change_indices_y = np.concatenate(([0], change_indices_y, [len(y) - 1]))
    
    # Getting the corresponding x and y values for both datasets
    x_steps_y = x[change_indices_y]
    y_steps = y[change_indices_y]

    x_steps_y = np.append(x_steps_y, x_steps_y[-1] + 1)  # Add 1 or an appropriate increment
    y_steps = np.append(y_steps, 0)
    

    # Creating the plot
    plt.figure(figsize=(12, 6), dpi=300)
    plt.step(x_steps_y, y_steps, where='post', label='# dereplicated covered genomes', color='blue')  # Specify color for y
    
    # Set title and labels
    plt.title(f'Histogram of the kmer coverage by their position | kmer size = {kmer_size}| {threshold} changes max')
    plt.xlabel('Kmer Start Position')
    plt.ylabel('Frequency')
    
    # Setting y-axis limits if needed
    # plt.ylim(-0.1, 1.1)
    
    # Add a legend
    plt.legend()

    plt.savefig(f'{save_path}/Histogram_of_the_kmer_coverage_{sample_name}_{kmer_size}.png', dpi=300)
    
    # Display the plot
    plt.show()
    
