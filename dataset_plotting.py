import pandas as pd
import matplotlib.pyplot as plt
import os

# Load the dataset
histology = pd.read_csv("./dataset/lesion_info.csv")

# Apply the histology mapping for adenoma vs. non-adenoma
histology_map = {
    "HP": "non-adenoma",
    "AD": "adenoma",
    "SSL": "non-adenoma",
    "TSA": "non-adenoma",
    "OTHER": "non-adenoma",
    "NO POLYP": "non-adenoma",
}
histology["histology_class"] = histology["histology_class"].replace(histology_map)

# Output directory for saving results
output_dir = "./histology_analysis"
os.makedirs(output_dir, exist_ok=True)

# 1. Overall Histogram: Number of Polyps per Histology Class
histology_counts = histology.groupby("histology_extended").size().reset_index(name="count")
histology_counts.to_csv(f"{output_dir}/histology_counts.csv", index=False)

plt.figure(figsize=(10, 6))
plt.bar(histology_counts["histology_extended"], histology_counts["count"], color='skyblue', edgecolor='black')
plt.title("Histogram of Number of Polyps per Histology Class")
plt.xlabel("Histology Class")
plt.ylabel("Number of Polyps")
plt.xticks(rotation=45, ha="right")
plt.grid(axis='y')
plt.tight_layout()
plt.savefig(f"{output_dir}/histology_histogram.png")
plt.close()

# 2. Overall Histogram: Adenoma vs Non-Adenoma
adenoma_counts = histology.groupby("histology_class").size().reset_index(name="count")

adenoma_counts.to_csv(f"{output_dir}/adenoma_counts.csv", index=False)

plt.figure(figsize=(10, 6))
plt.bar(adenoma_counts["histology_class"], adenoma_counts["count"], color=['salmon', 'lightgreen'], edgecolor='black')
plt.title("Histogram of Adenoma vs Non-Adenoma Polyps")
plt.ylabel("Number of Polyps")
plt.grid(axis='y')
plt.tight_layout()
plt.savefig(f"{output_dir}/adenoma_histogram.png")
plt.close()

# 3. Per Study Analysis: Histograms and Tables for Each Study
for video_name, group in histology.groupby("unique_video_name"):
    # Histology Class Distribution for this study
    study_histology_counts = group.groupby("histology_extended").size().reset_index(name="count")
    study_histology_counts.to_csv(f"{output_dir}/{video_name}_histology_counts.csv", index=False)

    plt.figure(figsize=(10, 6))
    plt.bar(study_histology_counts["histology_extended"], study_histology_counts["count"], color='skyblue', edgecolor='black')
    plt.title(f"Histology Distribution for {video_name}")
    plt.xlabel("Histology Class")
    plt.ylabel("Number of Polyps")
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{video_name}_histology_histogram.png")
    plt.close()

    # Adenoma vs Non-Adenoma Distribution for this study
    study_adenoma_counts = group.groupby("histology_class").size().reset_index(name="count")
    study_adenoma_counts.to_csv(f"{output_dir}/{video_name}_adenoma_counts.csv", index=False)

    plt.figure(figsize=(10, 6))
    plt.bar(study_adenoma_counts["histology_class"], study_adenoma_counts["count"], color=['salmon', 'lightgreen'], edgecolor='black')
    plt.title(f"Adenoma vs Non-Adenoma for {video_name}")
    plt.ylabel("Number of Polyps")
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{video_name}_adenoma_histogram.png")
    plt.close()

print("Histograms and tables have been saved to the output directory.")
