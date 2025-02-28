# REAL-Colon Dataset Analysis

## Overview

This directory contains scripts for analyzing the REAL-Colon dataset, focusing on lesion metadata and bounding box annotations. The analyses provide statistical insights and visualizations related to polyp sizes and frame-level annotation characteristics.

## Contents

### 1. polyp_histology.py

#### Description:
Computes and saves histograms of polyp histology in the REAL-Colon dataset.

#### Usage:

Update dataset_path = "./dataset/lesion_info.csv" with the correct path to the lesion metadata CSV file.
Update xml_folder = "./dataset/001-004_annotations" with the path to the folder containing XML annotation files.

Run the script using:

`python3 polyp_histology.py`

The histograms are saved to the output directory.

### 2. polyp_size.py

#### Description:
Computes and saves histograms of polyp sizes in the REAL-Colon dataset, categorized by histology classification.

#### Usage:

Update dataset_path = "./dataset/lesion_info.csv" with the correct path to the lesion metadata CSV file.

Run the script using:

`python3 polyp_size.py`

The histograms are saved to the output directory.

## Dependencies

Ensure the following Python libraries are installed before running the scripts:

`pip install matplotlib numpy pandas`

## Output

polyp_histology.py generates histograms of polyp histology and saves them in the output directory.

polyp_size.py generates histograms of polyp sizes and saves them in the output directory.
