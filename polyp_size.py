import os
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import glob

# Path to the directory where XML files are stored
xml_dir = "./dataset/001-004_annotations"
output_dir = "./output"
os.makedirs(output_dir, exist_ok=True)
histology = pd.read_csv("./dataset/lesion_info.csv")

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

# List of study IDs
histology["study_id"] = histology["unique_video_name"].str[:3]
study_ids = histology["study_id"].unique().tolist()

# Function to parse XML and extract polyp sizes
def extract_polyp_sizes_from_xml(xml_file, sizes):
    try:
        tree = ET.parse(xml_file)
        root = tree.getroot()
        size = root.find("./size")
        width, height = size.find("width"), size.find("height")
        
        if width is not None and height is not None:
            try:
                width = int(width.text.strip())
                height = int(height.text.strip())
                sizes.append(float(width * height))  # Store width * height directly
            except ValueError:
                print(f"Invalid size values in {xml_file}: width={width.text.strip()}, height={height.text.strip()}")
        else:
            print(f"Missing width/height in {xml_file}")
    except ET.ParseError:
        print(f"Error parsing {xml_file}")
    except Exception as e:
        print(f"Unexpected error with {xml_file}: {e}")
    
    return sizes

# Create a histogram of polyp sizes for each study
fig, axes = plt.subplots(2, 2, figsize=(12, 10), sharex=True)
axes = axes.flatten()

# Loop through studies
for i, study_id in enumerate(study_ids):
    ax = axes[i]
    processed_polyps = set()
    lesion_filtered = histology[histology["histology_class"].isin(["adenoma", "non-adenoma"])]

    # Filter for study_id and avoid redundant processing
    lesion_filtered_study = lesion_filtered[lesion_filtered["study_id"] == study_id]
    sizes_collected = []  # Collect all sizes for this study

    for _, row in lesion_filtered_study.iterrows():
        unique_object_id = row["unique_object_id"]
        if unique_object_id in processed_polyps:
            continue
        processed_polyps.add(unique_object_id)
        xml_files = sorted(glob.glob(os.path.join(xml_dir, f"{row['unique_video_name']}_*.xml")))
        
        if not xml_files:
            print(f"No XML files found for {row['unique_video_name']}")
            continue
    
        for xml_path in xml_files:
            # Extract polyp sizes from the XML file
            sizes = extract_polyp_sizes_from_xml(xml_path, sizes_collected)
            if sizes:  # Ensure sizes are added to the list
                sizes_collected.extend(sizes)

    # Only plot if there are sizes
    if sizes_collected:
        ax.hist(sizes_collected, bins=np.arange(0, max(sizes_collected) + 1, 1), edgecolor='black', color='skyblue', alpha=0.7)
    else:
        ax.set_title(f"No data for Study {study_id}")
        ax.set_xlabel("Size (mm)")
        ax.set_ylabel("Number of Polyps")

    # Set title and labels only after processing all XML files
    ax.set_title(f"Polyp Size Distribution for Study {study_id}")
    ax.set_xlabel("Size (mm)")
    ax.set_ylabel("Number of Polyps")
    ax.grid(True)

# Adjust layout to prevent overlap
fig.tight_layout(rect=[0, 0, 1, 0.96])

# Add a suptitle
fig.suptitle("Polyp Size Distribution Across Studies", fontsize=14)

# Save the plot
plt.savefig(f"{output_dir}/polyp_size_histograms.png")
plt.close()
