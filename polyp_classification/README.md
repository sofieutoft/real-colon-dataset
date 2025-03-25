# REAL-Colon Dataset Analysis

This directory contains scripts that provide statistical insights and visualizations related to the sizes and histology of polyps found in the REAL-Colon dataset.

## Contents

### 1. Polyp Histology

The script `python3 polyp_classification/polyp_histology.py` calculates and saves histograms depicting the distribution of polyp histologies in the REAL-Colon dataset. Before running the script, update the `dataset_path` variable to "./dataset/lesion_info.csv" to reflect the correct path to the lesion metadata CSV file. Similarly, update the `xml_folder` variable to "./dataset/001-004_annotations" for the path to the folder containing XML annotation files. The histograms are saved in the specified output directory.

### 2. Polyp Size

The script `python3 polyp_classification/polyp_size.py` generates and saves histograms of polyp sizes in the REAL-Colon dataset, categorized by histology classification. Ensure that the dataset paths are correctly set before execution to facilitate accurate data analysis.

## Description

overall_histology_histogram.png: Bar chart of polyps per histology class.

overall_avsna_histogram.png: Bar chart comparing adenoma vs. non-adenoma polyps.

combined_avsna_histograms.png: Grid of bar charts showing adenoma vs. non-adenoma distribution per study.

combined_histology_histograms.png: Grid of bar charts showing histology distribution per study.

adenoma_bbox_ratio_histogram.png: Histogram of bounding box ratios for adenomas.

non-adenoma_bbox_ratio_histogram.png: Histogram of bounding box ratios for non-adenomas.

polyp_size_histograms.png: Histogram showing the distribution of polyp sizes across the four cohorts.