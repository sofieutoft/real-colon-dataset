# REAL-Colon Dataset Analysis

This directory contains scripts that provide statistical insights and visualizations related to the sizes and histology of polyps found in the REAL-Colon dataset.

## Contents

### 1. Polyp Histology

The script `python3 polyp_histology.py` calculates and saves histograms depicting the distribution of polyp histologies in the REAL-Colon dataset. Before running the script, update the `dataset_path` variable to "./dataset/lesion_info.csv" to reflect the correct path to the lesion metadata CSV file. Similarly, update the `xml_folder` variable to "./dataset/001-004_annotations" for the path to the folder containing XML annotation files. The histograms are saved in the specified output directory.

### 2. Polyp Size

The script `python3 polyp_size.py` generates and saves histograms of polyp sizes in the REAL-Colon dataset, categorized by histology classification. Ensure that the dataset paths are correctly set before execution to facilitate accurate data analysis.