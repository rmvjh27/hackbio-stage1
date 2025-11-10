#a function for translating DNA to protein
def dna_to_protein(dna_sequence):
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'Stop', 'TAG': 'Stop', 'TGA': 'Stop'
    }
    protein = ""
    for i in range(0, len(dna_sequence) - 2, 3):  # Read in triplets
        codon = dna_sequence[i:i+3]
        amino_acid = codon_table.get(codon)
        if amino_acid == 'Stop':  # Stop codon
            break
        protein += amino_acid
    return protein

# a function for determining the Hamming distance between two strings
#making both strings same length
def hamming_distance(str1, str2):
  max_len = max(len(str1), len(str2))
  str1 = str1.ljust(max_len)
  str2 = str2.ljust(max_len)

  distance = 0
  for i in range(max_len):
    if str1[i] != str2[i]:
      distance += 1
  return distance

#usernames
slack_username = " Rama Jhowry"
twitter_handle = "@lord_chipo"

#calculate the hamming distance of Slack username and Twitter handle
distance = hamming_distance(slack_username, twitter_handle)

#print the output
print(f"The Hamming distance between'{slack_username}' and '{twitter_handle}' is {distance}.")

# Surprise task A
# Install and import the required libraries
!pip install pandas
!pip install seaborn
import pandas as pd
import seaborn as sns

# Import the task A part (a) dataset into the environment
data_source = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_top_deg_normalized_counts.csv" # provides URL to the data
df = pd.read_csv(data_source, index_col = 'Unnamed: 0') # loads the data into a pandas data frame

# visualize the data frame
df

# Create a clustered heatmap
sns.clustermap(df, cmap = 'Blues', figsize = (7,7))

# Import dataset for Task A part (b)
data_source2 = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/hbr_uhr_deg_chr22_with_significance.csv" # provides URL to the data
df_2 = pd.read_csv(data_source2, index_col = 'name') # loads the data into a pandas data frame

# visualize the data frame
df_2

# create a volcano plot from df_2
v_plot = sns.scatterplot(df_2, x = 'log2FoldChange', y = '-log10PAdj', hue = 'significance', palette = {'down':'orange', 'ns':'grey', 'up':'green'})

# add dashed vertical lines at log2FoldChange = ±1
import matplotlib.pyplot as plt
plt.axvline(x = 1, color = 'black', linestyle = '--')
plt.axvline(x = -1, color = 'black', linestyle = '--')

# Surprise task B
# Load the dataset from the uploaded file
from google.colab import files
uploaded = files.upload()
df = pd.read_csv('Breast Cancer Wisconsin Dataset (1).txt')

print(df.head(5))

print(df.shape)

"""# **Scatter Plot**"""
plt.figure(figsize=(8, 6))
sns.scatterplot(data=df,
                x="radius_mean",
                y="texture_mean",
                hue="diagnosis")

"""## **Heatmap**"""

# Select the specified columns
selected_cols = [
    'radius_mean', 'texture_mean', 'perimeter_mean',
    'area_mean', 'smoothness_mean', 'compactness_mean'
]
df_subset = df[selected_cols]

# Compute the correlation matrix
corr_matrix = df_subset.corr()

# Create a heatmap of the correlation matrix
plt.figure(figsize=(10, 8))
sns.heatmap(corr_matrix, annot=True, cmap='Blues', fmt=".1f")

# Breast Cancer Diagnostic Data Visualization
# Import dataset from the provided URL

url = "https://raw.githubusercontent.com/HackBio-Internship/2025_project_collection/refs/heads/main/Python/Dataset/data-3.csv"

# Read the dataset directly from GitHub
data = pd.read_csv(url)

# Display the first 5 rows to confirm successful import
print(" Dataset successfully loaded!\n")
print(data.head())

# Basic info about the dataset

print(" \n Dataset Information:")
print(data.info())

# Check for missing values
print("\n Missing Values in Each Column:")
print(data.isnull().sum())

# Scatter Plot (Smoothness vs Compactness)
# Purpose:
# Visualize relationship between 'smoothness_mean' and 'compactness_mean'
# Color the data points based on the diagnosis (M = Malignant, B = Benign)

# Set a consistent Seaborn style
sns.set(style="whitegrid")

# Define figure size
plt.figure(figsize=(8,6))

# Create scatter plot
sns.scatterplot(
    data=data,
    x="compactness_mean",        # X-axis variable
    y="smoothness_mean",         # Y-axis variable
    hue="diagnosis",             # Color points by diagnosis type
    palette="Set1",              # Distinct colors for M and B
    s=70,                        # Size of scatter points
    edgecolor="black"            # Black border for clarity
)

# Add title and axis labels
plt.title("Scatter Plot: Smoothness vs Compactness", fontsize=14)
plt.xlabel("Compactness (Mean)", fontsize=12)
plt.ylabel("Smoothness (Mean)", fontsize=12)

# Add gridlines for better readability
plt.grid(True, linestyle="--", alpha=0.6)

# Show the plot
plt.show()

# Density Plot (Area Distribution)
# Purpose:
# Plot the distribution of 'area_mean' for both Malignant and Benign cases
# using smooth Kernel Density Estimates (KDE)

# Define figure size
plt.figure(figsize=(8,6))

# KDE for Malignant (M)
sns.kdeplot(
    data=data[data["diagnosis"] == "M"]["area_mean"],  # Subset only M cases
    label="Malignant",         # Label for the legend
    fill=True,                 # Fill area under curve
    alpha=0.5,                 # Transparency for overlap clarity
    linewidth=2
)

# KDE for Benign (B)
sns.kdeplot(
    data=data[data["diagnosis"] == "B"]["area_mean"],  # Subset only B cases
    label="Benign",            # Label for the legend
    fill=True,
    alpha=0.5,
    linewidth=2
)

# Add title, axis labels, and legend
plt.title("Density Plot of Area Mean by Diagnosis", fontsize=14)
plt.xlabel("Area (Mean)", fontsize=12)
plt.ylabel("Density", fontsize=12)
plt.legend(title="Diagnosis")

# Add gridlines for readability
plt.grid(True, linestyle="--", alpha=0.6)

# Display the plot
plt.show()

