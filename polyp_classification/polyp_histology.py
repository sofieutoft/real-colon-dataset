"""
Compute and save histograms of the histology of the polyps in the REAL-Colon dataset.

Usage:
    - Update dataset_path = "./dataset/lesion_info.csv" with the path to the lesion metadata CSV file.
    - Update xml_folder = "./dataset/001-004_annotations" with the path to the folder containing XML annotation files.
    - python3 polyp_histology.py

"""

import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET

# ---------------------- Load and Preprocess Dataset ----------------------

def load_histology_data(filepath):
    histology = pd.read_csv(filepath)
    histology_map = {
        "HP": "non-adenoma",
        "AD": "adenoma",
        "SSL": "non-adenoma",
        "TSA": "non-adenoma",
        "OTHER": "non-adenoma",
        "NO POLYP": "non-adenoma",
    }
    histology["histology_class"] = histology["histology_class"].replace(histology_map)
    histology["study_id"] = histology["unique_video_name"].str[:3]
    return histology

# ---------------------- Helper Function to Save Histograms ----------------------

def save_histogram(data, x_col, title, xlabel, ylabel, filename, colors=None):
    plt.figure(figsize=(10, 6))
    plt.bar(data[x_col], data["count"], color=colors, edgecolor='black')
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()

# ---------------------- Generate Histograms ----------------------

def generate_histograms(histology, output_dir):
    os.makedirs(output_dir, exist_ok=True)
    
    histology_counts = histology.groupby("histology_extended").size().reset_index(name="count")
    histology_counts.to_csv(f"{output_dir}/overall_histology_counts.csv", index=False)
    save_histogram(histology_counts, "histology_extended", "Histogram of Polyps per Histology Class", "Histology Class", "Number of Polyps", f"{output_dir}/overall_histology_histogram.png", colors='skyblue')
    
    adenoma_counts = histology.groupby("histology_class").size().reset_index(name="count")
    adenoma_counts.to_csv(f"{output_dir}/overall_adenoma_counts.csv", index=False)
    save_histogram(adenoma_counts, "histology_class", "Adenoma vs Non-Adenoma Polyps", "Histology Class", "Number of Polyps", f"{output_dir}/overall_adenoma_histogram.png", colors=['salmon', 'lightgreen'])

    # Generate per-study adenoma histograms in a 2x2 grid
    study_ids = histology["study_id"].unique()
    n_studies = len(study_ids)
    cols = 2
    rows = (n_studies + 1) // cols

    fig, axes = plt.subplots(rows, cols, figsize=(10, 6 * rows), sharey=True)
    axes = axes.flatten()

    for ax, study_id in zip(axes, study_ids):
        group = histology[histology["study_id"] == study_id]
        study_adenoma_counts = group.groupby("histology_class").size().reset_index(name="count")
        ax.bar(study_adenoma_counts["histology_class"], study_adenoma_counts["count"], alpha=0.7, edgecolor='black')
        ax.set_title(f"Cohort {study_id}")
        ax.set_xlabel("Histology Class")
        ax.set_xticklabels(study_adenoma_counts["histology_class"], rotation=45, ha="right")
        ax.grid(axis='y')

    # Hide unused subplots
    for ax in axes[len(study_ids):]:
        ax.axis('off')

    fig.suptitle("Adenoma vs Non-Adenoma Polyps Across Cohorts")
    fig.tight_layout()
    plt.savefig(f"{output_dir}/combined_adenoma_histograms.png")
    plt.close()

    # Generate per-study histology histograms
    fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=False)

    # Flatten axes array for easy iteration
    axes = axes.flatten()

    for i, study_id in enumerate(study_ids):
        ax = axes[i]  # Get the corresponding subplot
        
        # Filter data for the current study
        group = histology[histology["study_id"] == study_id]
        study_histology_counts = group.groupby("histology_extended").size().reset_index(name="count")
        
        # Create bar plot
        ax.bar(study_histology_counts["histology_extended"], study_histology_counts["count"], edgecolor='black')
        
        ax.set_title(f"Histology Distribution for Cohort {study_id}")
        ax.set_ylabel("Number of Polyps")
        ax.grid(axis='y')

        # Fix x-tick labels
        ax.set_xticks(study_histology_counts.index)  # Use actual index positions
        ax.set_xticklabels(study_histology_counts["histology_extended"], rotation=45, ha="right")

    # Remove any unused subplots (if fewer than 4 studies exist)
    for j in range(i + 1, len(axes)):
        fig.delaxes(axes[j])

    fig.suptitle("Histology Distribution Across Cohorts", fontsize=14)
    fig.tight_layout(rect=[0, 0, 1, 0.96])  # Leave space for the title
    plt.savefig(f"{output_dir}/combined_histology_histograms.png")
    plt.close()


# ---------------------- Extract Bounding Box Ratios ----------------------

def get_bbox_ratio(xml_path):
    try:
        tree = ET.parse(xml_path)
        root = tree.getroot()
        size = root.find("./size")
        width, height = int(size.find("width").text), int(size.find("height").text)
        image_diag = (width ** 2 + height ** 2) ** 0.5
        
        obj = root.find("./object")
        bbox = obj.find("bndbox")
        xmin, xmax, ymin, ymax = [int(bbox.find(tag).text) for tag in ["xmin", "xmax", "ymin", "ymax"]]
        bbox_diag = ((xmax - xmin) ** 2 + (ymax - ymin) ** 2) ** 0.5
        return bbox_diag / image_diag
    except:
        return None

# ---------------------- Process Annotations ----------------------

def process_annotations(histology, xml_folder):
    bbox_ratios = []
    processed_polyps = set()
    lesion_filtered = histology[histology["histology_class"].isin(["adenoma", "non-adenoma"])]
    
    for _, row in lesion_filtered.iterrows():
        unique_object_id = row["unique_object_id"]
        if unique_object_id in processed_polyps:
            continue
        processed_polyps.add(unique_object_id)
        
        xml_files = sorted(glob.glob(os.path.join(xml_folder, f"{row['unique_video_name']}_*.xml")))
        if not xml_files:
            continue
        
        for xml_path in xml_files:
            bbox_ratio = get_bbox_ratio(xml_path)
            if bbox_ratio is not None:
                bbox_ratios.append({
                    "unique_object_id": unique_object_id,
                    "unique_video_name": row["unique_video_name"],
                    "size_mm": row["size [mm]"],
                    "site": row["site"],
                    "histology_extended": row["histology_extended"],
                    "histology_class": row["histology_class"],
                    "bbox_ratio": bbox_ratio
                })
    return pd.DataFrame(bbox_ratios)

# ---------------------- Generate Bounding Box Histograms ----------------------

def generate_bbox_histograms(bbox_df, output_dir):
    for hist_class in ["adenoma", "non-adenoma"]:
        subset = bbox_df[bbox_df["histology_class"] == hist_class]
        if subset.empty:
            print(f"Warning: No data for {hist_class}, skipping histogram.")
            continue
        
        plt.figure(figsize=(10, 6))
        plt.hist(subset["bbox_ratio"], bins=20, color='blue' if hist_class == "adenoma" else 'green', edgecolor='black')
        plt.title(f"Histogram of BBox Ratio for {hist_class.capitalize()} Polyps")
        plt.xlabel("Ratio of BBox Diagonal to Image Diagonal")
        plt.ylabel("Number of Frames")
        plt.grid(axis='y')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{hist_class}_bbox_ratio_histogram.png")
        plt.close()

# ---------------------- Main Execution ----------------------

def main():
    dataset_path = "./dataset/lesion_info.csv"
    xml_folder = "./dataset/001-004_annotations"
    output_dir = "./output"
    
    histology = load_histology_data(dataset_path)
    generate_histograms(histology, output_dir)
    
    bbox_df = process_annotations(histology, xml_folder)
    bbox_df.to_csv(f"{output_dir}/bbox_ratios.csv", index=False)
    
    generate_bbox_histograms(bbox_df, output_dir)
    print("Histograms and tables have been saved to the output directory.")

if __name__ == "__main__":
    main()
