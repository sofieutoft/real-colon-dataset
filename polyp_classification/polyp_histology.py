"""
Compute and save histograms of the histology of the polyps in the REAL-Colon dataset.

Usage:
    - Update dataset_path = "./dataset/lesion_info.csv" with the path to the lesion metadata CSV file.
    - Update xml_folder = "./dataset/001-004_annotations" with the path to the folder containing XML annotation files.
    - python3 polyp_classification/polyp_histology.py

"""
import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import xml.etree.ElementTree as ET
import concurrent.futures
import argparse

def load_histology_data(filepath):
    """
    Load the histology data from the given CSV file and map the histology classes to adenoma and non-adenoma.

    Args:
        filepath (str): The path to the CSV file containing the histology data.
    
    Returns:
        pd.DataFrame: The histology data with the histology classes mapped to adenoma and non-adenoma.
    """
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

def save_histogram(data, x_col, title, xlabel, ylabel, filename, colors=None):
    """
    Save a histogram of the given data to the specified file.
    
    Args:
        data (pd.DataFrame): The data to plot.
        x_col (str): The column to use for the x-axis.
        title (str): The title of the plot.
        xlabel (str): The label for the x-axis.
        ylabel (str): The label for the y-axis.
        filename (str): The path to save the plot.
        colors (list): The colors to use for the bars.
    
    Generates:
        A histogram of the given data saved to the specified file.
    """
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

def generate_histograms(histology, output_dir, save_csv=True):
    """
    Generate and save histograms based on the histology data.

    Args:
        histology (pd.DataFrame): The histology data to analyze.
        output_dir (str): The directory to save the output files.
        save_csv (bool): Whether to save the histograms to CSV files.
    
    Generates:
        Histograms and optionally CSV files of the histology data.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    histology_counts = histology.groupby("histology_extended").size().reset_index(name="count")
    if save_csv:
        histology_counts.to_csv(f"{output_dir}/overall_histology_counts.csv", index=False)
    save_histogram(histology_counts, "histology_extended", "Histogram of Polyps per Histology Class", "Histology Class", "Number of Polyps", f"{output_dir}/overall_histology_histogram.png", colors='skyblue')
    
    adenoma_counts = histology.groupby("histology_class").size().reset_index(name="count")
    if save_csv:
        adenoma_counts.to_csv(f"{output_dir}/overall_adenoma_counts.csv", index=False)
    save_histogram(adenoma_counts, "histology_class", "Adenoma vs Non-Adenoma Polyps", "Histology Class", "Number of Polyps", f"{output_dir}/overall_avsna_histogram.png", colors=['salmon', 'lightgreen'])

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
    plt.savefig(f"{output_dir}/combined_avsna_histograms.png")
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


def get_bbox_ratio(xml_path):
    """
    Calculate the ratio of the bounding box diagonal to the image diagonal.

    Args:
        xml_path (str): The path to the XML file containing the bounding box and image size information.
    
    Returns:
        float: The ratio of the bounding box diagonal to the image diagonal, or None if an error occurs.
    """
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


def process_annotations(histology, xml_folder):
    """
    Process the XML annotation files to calculate bounding box ratios for each polyp.

    Args:
        histology (pd.DataFrame): The histology data containing polyp information.
        xml_folder (str): The path to the folder containing XML annotation files.
    
    Returns:
        pd.DataFrame: A DataFrame containing the bounding box ratios for each polyp.
    """
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


def generate_bbox_histograms(bbox_df, output_dir):
    """
    Generates and saves histograms of bounding box ratios for different histology classes.
    Args:
        bbox_df (pd.DataFrame): DataFrame containing bounding box data with a column 'histology_class' 
                                indicating the class of the polyp and a column 'bbox_ratio' indicating 
                                the ratio of the bounding box diagonal to the image diagonal.
        output_dir (str): Directory where the histogram images will be saved.
    Raises:
        Warning: If there is no data for a specific histology class, a warning message is printed and 
                 the histogram for that class is skipped.
    Returns:
        None
    """
    for hist_class in ["adenoma", "non-adenoma"]:

        subset = bbox_df[bbox_df["histology_class"] == hist_class]
        if subset.empty:
            print(f"Warning: No data for {hist_class}, skipping histogram.")
            continue

        plt.figure(figsize=(10, 6))
        plt.hist(subset["bbox_ratio"], bins=20, color='blue' if hist_class == "adenoma" else 'green', edgecolor='black')
        plt.title(f"Histogram of BBox Ratio for {hist_class.capitalize()} Polyps")
        plt.xlabel("Ratio of BBox Diagonal to Image Diagonal")
        plt.ylabel("Number aof Frames")
        plt.grid(axis='y')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/{hist_class}_bbox_ratio_histogram.png")
        plt.close()


def main(save_csv):
    """
    Processes histology data and generates histograms.

    Args:
        save_csv (bool): If True, saves bounding box ratios as a CSV file.

    This function:
    - Loads histology data from the dataset path.
    - Extracts bounding box data from annotation folders.
    - Generates and saves histograms of bounding box data.
    - Saves extracted data as a CSV file if `save_csv` is True.
    """
    dataset_path = "/ssd/storage/shared/colonscopy/public_datasets/real-colon_dataset_released_v20230228/"
    output_dir = "./stats"
    os.makedirs(output_dir, exist_ok=True)

    histology = load_histology_data(os.path.join(dataset_path, "lesion_info.csv"))
    generate_histograms(histology, output_dir, save_csv)

    annotation_folders = [
        os.path.join(dataset_path, folder)
        for folder in os.listdir(dataset_path)
        if folder.endswith("_annotations") and os.path.isdir(os.path.join(dataset_path, folder))
    ]

    # Helper to wrap the process_annotations call with histology
    def process_folder(folder_path):
        df = process_annotations(histology, folder_path)
        return df if not df.empty else None

    with concurrent.futures.ProcessPoolExecutor(max_workers=16) as executor:
        results = list(executor.map(process_folder, annotation_folders))

    bbox_dfs = [df for df in results if df is not None]

    if bbox_dfs:
        bbox_df = pd.concat(bbox_dfs, ignore_index=True)
        if save_csv:
            bbox_df.to_csv(f"{output_dir}/bbox_ratios.csv", index=False)
        generate_bbox_histograms(bbox_df, output_dir)
    else:
        print("No annotation data found in any _annotations folder.")

    print("Histograms and tables have been saved to the output directory.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Save to CSV")
    parser.add_argument("--CSV", default=True, help="Save histology analysis to CSV files.")
    args = parser.parse_args()
    main(args.CSV)
