import pandas as pd
import matplotlib.pyplot as plt
import os
import xml.etree.ElementTree as ET
import glob


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

# Extract study ID from unique_video_name
histology["study_id"] = histology["unique_video_name"].str[:3]

# Output directory for saving results
output_dir = "./histology_analysis"
os.makedirs(output_dir, exist_ok=True)

# 1. Overall Histogram: Number of Polyps per Histology Class
histology_counts = histology.groupby("histology_extended").size().reset_index(name="count")
histology_counts.to_csv(f"{output_dir}/overall_histology_counts.csv", index=False)

plt.figure(figsize=(10, 6))
plt.bar(histology_counts["histology_extended"], histology_counts["count"], color='skyblue', edgecolor='black')
plt.title("Histogram of Number of Polyps per Histology Class")
plt.xlabel("Histology Class")
plt.ylabel("Number of Polyps")
plt.xticks(rotation=45, ha="right")
plt.grid(axis='y')
plt.tight_layout()
plt.savefig(f"{output_dir}/overall_histology_histogram.png")
plt.close()

# 2. Overall Histogram: Adenoma vs Non-Adenoma
adenoma_counts = histology.groupby("histology_class").size().reset_index(name="count")
adenoma_counts.to_csv(f"{output_dir}/overall_adenoma_counts.csv", index=False)

plt.figure(figsize=(10, 6))
plt.bar(adenoma_counts["histology_class"], adenoma_counts["count"], color=['salmon', 'lightgreen'], edgecolor='black')
plt.title("Histogram of Adenoma vs Non-Adenoma Polyps")
plt.ylabel("Number of Polyps")
plt.grid(axis='y')
plt.tight_layout()
plt.savefig(f"{output_dir}/overall_adenoma_histogram.png")
plt.close()

# 3. Per Study Analysis: Histograms and Tables for Each Study
for study_id, group in histology.groupby("study_id"):
    # Histology Class Distribution for this study
    study_histology_counts = group.groupby("histology_extended").size().reset_index(name="count")
    study_histology_counts.to_csv(f"{output_dir}/{study_id}_histology_counts.csv", index=False)

    plt.figure(figsize=(10, 6))
    plt.bar(study_histology_counts["histology_extended"], study_histology_counts["count"], color='skyblue', edgecolor='black')
    plt.title(f"Histology Distribution for Study {study_id}")
    plt.xlabel("Histology Class")
    plt.ylabel("Number of Polyps")
    plt.xticks(rotation=45, ha="right")
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{study_id}_histology_histogram.png")
    plt.close()

    # Adenoma vs Non-Adenoma Distribution for this study
    study_adenoma_counts = group.groupby("histology_class").size().reset_index(name="count")
    study_adenoma_counts.to_csv(f"{output_dir}/{study_id}_adenoma_counts.csv", index=False)

    plt.figure(figsize=(10, 6))
    plt.bar(study_adenoma_counts["histology_class"], study_adenoma_counts["count"], color=['salmon', 'lightgreen'], edgecolor='black')
    plt.title(f"Adenoma vs Non-Adenoma for Study {study_id}")
    plt.ylabel("Number of Polyps")
    plt.grid(axis='y')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/{study_id}_adenoma_histogram.png")
    plt.close()

def get_bbox_ratio(xml_path):
    """Extracts the bounding box ratio from an XML file."""
    tree = ET.parse(xml_path)
    root = tree.getroot()
    
    size = root.find("./size")
    if size is None:
        return None
    
    width = size.find("width")
    height = size.find("height")
    if width is None or height is None:
        return None
    
    image_diag = ((int(width.text) ** 2 + int(height.text) ** 2) ** 0.5)
    obj = root.find("./object")
    if obj is None:
        return None
    
    bbox = obj.find("bndbox")
    if bbox is None:
        return None
    
    coords = [bbox.find(tag) for tag in ["xmin", "xmax", "ymin", "ymax"]]
    if any(coord is None for coord in coords):
        return None
    
    xmin, xmax, ymin, ymax = [int(coord.text) for coord in coords]
    bbox_diag = (((xmax - xmin) ** 2 + (ymax - ymin) ** 2) ** 0.5)
    return bbox_diag / image_diag

def process_annotations(lesion_info, xml_folder):
    """Processes bounding box ratios while ensuring each polyp (`unique_object_id`) is processed only once."""
    bbox_ratios = []
    processed_polyps = set()

    # Filter to only adenoma and non-adenoma cases
    lesion_filtered = lesion_info[lesion_info["histology_class"].isin(["adenoma", "non-adenoma"])]

    for _, row in lesion_filtered.iterrows():
        unique_object_id = row["unique_object_id"]
        if unique_object_id in processed_polyps:
            continue  # Skip duplicate polyps
        
        processed_polyps.add(unique_object_id)
        
        xml_pattern = os.path.join(xml_folder, f"{row['unique_video_name']}_*.xml")
        xml_files = sorted(glob.glob(xml_pattern))

        if not xml_files:
            continue  # Skip if no annotation files found
        
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


bbox_df = process_annotations(histology, "./dataset/001-004_annotations")
bbox_df.to_csv(f"{output_dir}/bbox_ratios.csv", index=False)

# Plot and save histograms for each histology class
for hist_class in ["adenoma", "non-adenoma"]:
    subset = bbox_df[bbox_df["histology_class"] == hist_class]

    if subset.empty:
        print(f"Warning: No data for {hist_class}, skipping histogram.")
        continue

    plt.figure(figsize=(10, 6))
    plt.hist(subset["bbox_ratio"], bins=20, color='blue' if hist_class == "adenoma" else 'green', edgecolor='black')
    plt.title(f"Histogram of BBox Ratio for {hist_class} Polyps")
    plt.xlabel("Ratio of BBox Diagonal to Image Diagonal")
    plt.ylabel("Number of Frames")
    plt.grid(axis='y')
    plt.tight_layout()

    output_path = f"{output_dir}/{hist_class}_bbox_ratio_histogram.png"
    plt.savefig(output_path)
    plt.close()

print("Histograms and tables have been saved to the output directory.")
