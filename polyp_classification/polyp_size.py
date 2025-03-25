"""
Compute and save histograms of the size of the polyps in the REAL-Colon dataset.

Usage:
    - Update dataset_path = "./dataset/lesion_info.csv" with the path to the lesion metadata CSV file.
    - python3 polyp_classification/polyp_size.py

"""

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

def main():
    # Define output directory
    output_dir = "./stats"
    os.makedirs(output_dir, exist_ok=True)

    # Load dataset
    histology = pd.read_csv("/ssd/storage/shared/colonscopy/public_datasets/real-colon_dataset_released_v20230228/lesion_info.csv")

    # Define histology mapping for adenoma vs. non-adenoma
    histology_map = {
        "HP": "non-adenoma",
        "AD": "adenoma",
        "SSL": "non-adenoma",
        "TSA": "non-adenoma",
        "OTHER": "non-adenoma",
        "NO POLYP": "non-adenoma",
    }
    histology["histology_class"] = histology["histology_class"].replace(histology_map)

    # Extract study IDs
    histology["study_id"] = histology["unique_video_name"].str[:3]
    study_ids = histology["study_id"].unique().tolist()

    # Create histograms for polyp sizes
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
    axes = axes.flatten()

    for i, study_id in enumerate(study_ids):
        ax = axes[i]
        processed_polyps = set()
        lesion_filtered = histology[histology["histology_class"].isin(["adenoma", "non-adenoma"])]
        lesion_filtered_study = lesion_filtered[lesion_filtered["study_id"] == study_id]
        sizes_collected = []

        for _, row in lesion_filtered_study.iterrows():
            unique_object_id = row["unique_object_id"]
            if unique_object_id in processed_polyps:
                continue
            processed_polyps.add(unique_object_id)
            size = row["size [mm]"]
            if pd.notna(size):
                sizes_collected.append(float(size))

        # Plot histogram if data exists
        if sizes_collected:
            ax.hist(sizes_collected, bins=np.arange(0, max(sizes_collected) + 1, 1), edgecolor='black', color='skyblue', alpha=0.7)
        else:
            ax.set_title(f"No data for Study {study_id}")

        ax.set_title(f"Polyp Size Distribution for Cohort {study_id}")
        ax.set_xlabel("Size (mm)")
        ax.set_ylabel("Number of Polyps")
        ax.grid(True)

    # Adjust layout
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    fig.suptitle("Polyp Size Distribution Across Cohorts", fontsize=14)

    # Save plot
    plt.savefig(f"{output_dir}/polyp_size_histograms.png")
    plt.close()
    print("Histograms have been saved to the output directory.")

if __name__ == "__main__":
    main()