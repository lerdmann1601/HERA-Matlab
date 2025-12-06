"""
HERA_EXAMPLE_2_WORKFLOW - Example workflow for HERA analysis with pruning and stability checking.

Description:
  This script demonstrates a complete workflow for running HERA on a series of datasets 
  (e.g., decreasing contrast/alpha levels). It performs the following steps:
  1. Sets up the necessary directory structure.
  2. Iterates through data folders (sorted by alpha).
  3. Runs the HERA analysis (via compiled runtime or local script).
  4. Analyzes the results, specifically calculating Confidence Intervals (CIs) from bootstrap ranks.
  5. Prunes the results to retain only the transition point and stable regions where methods overlap.
  6. Generates a stability curve plot.

Workflow:
  1. Initialization: Define paths and configuration.
  2. Setup: Ensure data and results directories exist.
  3. Data Processing: 
     - Detect input folders.
     - Configure HERA for each dataset.
     - Execute HERA.
     - specific analysis of JSON output.
  4. Pruning & Plotting:
     - Identify overlap regions between methods (e.g., Method G and Method B).
     - Keep relevant datasets (overlap + transition point).
     - Remove redundant stable datasets.
     - Plot the optimization curve.

Usage:
  Run this script directly from the terminal or IDE:
  $ python3 hera_example_2_workflow.py

Author: Lukas von Erdmannsdorff
"""

import os
import sys
import shutil
import json
import subprocess
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

# --- Configuration ---
BASE_DIR = Path(__file__).parent.resolve()
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
HERA_DEPLOY_DIR = BASE_DIR.parent / "release" / "macos"
HERA_EXECUTABLE_NAME = "run_HERA_Runtime.sh"
HERA_MCR_PATH = "/Applications/MATLAB_R2025b.app"

METRIC_NAMES = ["OC", "CNR", "SNR"]
RANKING_MODE = "M1_M2_M3"

# Manual Bootstrap Settings
MANUAL_BOOTSTRAP = {
    "manual_B_thr": 2000,
    "manual_B_ci": 10000,
    "manual_B_rank": 500
}

class HERAWorkflow:
    def __init__(self):
        self.data_dir = DATA_DIR
        self.results_dir = RESULTS_DIR
        
    def setup_directories(self):
        """
        1. Ensures the directory structure exists.
        Checks for data and result directories and creates them if missing.
        """
        if not self.data_dir.exists():
            print(f"[Info] Creating Data Directory: {self.data_dir}")
            self.data_dir.mkdir(parents=True)
            
        if not self.results_dir.exists():
            print(f"[Info] Creating Results Directory: {self.results_dir}")
            self.results_dir.mkdir(parents=True)
            
    # migrate_legacy_data removed per user request


            
    def create_hera_config(self, input_path, output_path):
        """
        3. Creates the JSON configuration structure for HERA.
        Constructs the dictionary expected by the HERA runtime.
        """
        config = {
            "userInput": {
                # --- Mandatory ---
                "folderPath": str(input_path),
                "metric_names": METRIC_NAMES,
                "output_dir": str(output_path),
                
                # --- Optional ---
                "fileType": ".csv",
                "ranking_mode": RANKING_MODE,
                "reproducible": True,
                "create_reports": False,  # Disabled for speed as requested
                "ci_level": 0.95,
                
                # --- Manual Bootstrap ---
                "manual_B_thr": MANUAL_BOOTSTRAP["manual_B_thr"],
                "manual_B_ci": MANUAL_BOOTSTRAP["manual_B_ci"],
                "manual_B_rank": MANUAL_BOOTSTRAP["manual_B_rank"]
            }
        }
        return config

    def run_hera(self, config_path):
        """
        4. Executes HERA using the compiled runtime/launcher.
        Locates the executable and runs it with the provided configuration.
        """
        hera_cmd = HERA_DEPLOY_DIR / HERA_EXECUTABLE_NAME
        
        if not hera_cmd.exists():
             # Fallback check
            local_cmd = BASE_DIR.parent / "HERA_Launcher.command"
            if local_cmd.exists():
                 hera_cmd = local_cmd
            else:
                 print(f"[Error] HERA executable not found at {hera_cmd}")
                 return False

        print(f" -> Running HERA with config: {config_path}")
        
        cmd = [str(hera_cmd), HERA_MCR_PATH, "configFile", str(config_path)]
        
        try:
            result = subprocess.run(cmd, capture_output=True, text=True)
            if result.returncode != 0:
                print("[Error] HERA Run Failed!")
                print(result.stderr)
                return False
            return True
        except Exception as e:
            print(f"[Error] Execution Exception: {e}")
            return False

    def get_latest_result_json(self, result_folder):
        """
        Helper: Finds the latest JSON output in the results folder.
        """
        if not result_folder.exists():
            return None
        
        # Sort by creation time (implicitly by name if timestamped)
        timestamp_dirs = sorted([d for d in result_folder.iterdir() if d.is_dir()])
        if not timestamp_dirs:
            return None
            
        latest = timestamp_dirs[-1]
        output_dir = latest / "Output"
        if not output_dir.exists():
            return None
            
        jsons = list(output_dir.glob("*.json"))
        if not jsons:
            return None
            
        return jsons[0]

    def analyze_results(self, json_path):
        """
        5. Analysis: Reads ranking info and calculates CIs from bootstrap ranks.
        Specific logic to extract bootstrap distributions and compute percentiles.
        """
        with open(json_path, 'r') as f:
            data = json.load(f)
            
        results = data.get('results', {})
        dataset_names = data.get('dataset_names', [])
        
        final_ranks = results.get('final_rank', [])
        
        # KEY CHANGE: Calculate CIs from bootstrap distribution
        # Matlab 'final_bootstrap_ranks' is [num_datasets x B]
        # JSON export usually preserves this as a list of lists (rows)
        bootstrap_ranks = results.get('final_bootstrap_ranks', [])
        
        ci_lower = []
        ci_upper = []
        
        if bootstrap_ranks and len(bootstrap_ranks) == len(dataset_names):
            try:
                # Convert to numpy array for easy percentile calc
                # Expected shape: (Num_Datasets, B)
                bs_data = np.array(bootstrap_ranks)
                
                # Calculate 95% CI (2.5th and 97.5th percentiles)
                # axis=1 means along the B dimension (columns)
                ci_lower = np.percentile(bs_data, 2.5, axis=1)
                ci_upper = np.percentile(bs_data, 97.5, axis=1)
                print(f" -> Calculated CIs from {bs_data.shape[1]} bootstrap samples.")
            except Exception as e:
                print(f"[Warning] Error calculating CIs from bootstrap ranks: {e}")
                ci_lower = final_ranks
                ci_upper = final_ranks
        else:
             print(f"[Warning] Bootstrap ranks missing or shape mismatch in {json_path.name}")
             ci_lower = final_ranks
             ci_upper = final_ranks

        rank_map = {}
        for i, name in enumerate(dataset_names):
            rank_map[name] = {
                "rank": final_ranks[i],
                "ci_min": ci_lower[i],
                "ci_max": ci_upper[i]
            }
            
        return rank_map

    def process_data(self):
        """
        Main Loop: Iterates through data folders, runs HERA, and collects results.
        """
        self.setup_directories()


        
        # Find all Alpha_XX folders in python_experiments/data
        alpha_folders = [d for d in self.data_dir.iterdir() if d.is_dir() and d.name.startswith("Alpha_")]
        
        # Sort by alpha value descending (High Alpha/Contrast -> Low Alpha/Noise)
        def parse_alpha(p):
            try:
                return int(p.name.split("_")[-1])
            except:
                return -1
                
        alpha_folders.sort(key=parse_alpha, reverse=True)
        
        if not alpha_folders:
            print("[Info] No Alpha_XX data folders found to process.")
            return

        aggregated_results = []

        print("\n--- Starting Processing Loop ---")
        for data_folder in alpha_folders:
            alpha_str = data_folder.name # Alpha_XX
            alpha_val = parse_alpha(data_folder)
            
            result_target_dir = self.results_dir / alpha_str
            if not result_target_dir.exists():
                result_target_dir.mkdir(parents=True)
                
            config_path = result_target_dir / "config.json"
            
            # Create Config
            config_dict = self.create_hera_config(data_folder, result_target_dir)
            with open(config_path, 'w') as f:
                json.dump(config_dict, f, indent=4)
                
            # Run HERA
            # Strategy: Always run to guarantee consistent state with CIs.
            success = self.run_hera(config_path)
            # success = True # Skip for fast verification using cached results
            if not success:
                print(f"[Skip] Failed to process {alpha_str}")
                continue
                
            # Analyze Result
            json_path = self.get_latest_result_json(result_target_dir)
            if not json_path:
                print(f"[Skip] No output JSON found for {alpha_str}")
                continue
                
            stats = self.analyze_results(json_path)
            
            res_entry = {
                "Alpha": alpha_val / 100.0,
                "Folder": result_target_dir,
                "Stats": stats
            }
            aggregated_results.append(res_entry)
            
        self.prune_and_plot(aggregated_results)

    def prune_and_plot(self, results):
        """
        6. Pruning and Plotting Logic.
        
        Pruning Rules:
          1. Keep all results where Method G and Method B CIs OVERLAP.
          2. Keep the method immediately BEFORE the overlap starts (Transition point).
          3. Delete all others (Stable results).
        """
        if not results:
            return
            
        print("\n--- Pruning and Plotting ---")
        
        # Sort results by Alpha Descending (High -> Low)
        results.sort(key=lambda x: x["Alpha"], reverse=True)
        
        # Identify Overlaps
        to_keep_indices = set()
        
        first_overlap_index = -1
        
        for i, res in enumerate(results):
            g_stats = res["Stats"].get("Method G", {})
            b_stats = res["Stats"].get("Method B", {})
            
            # Ensure both methods exist and have CI data
            if not (g_stats and b_stats and "ci_min" in g_stats and "ci_max" in g_stats and "ci_min" in b_stats and "ci_max" in b_stats):
                print(f" [Warning] Skipping overlap check for Alpha {res['Alpha']} due to missing Method G or B stats.")
                continue

            # Check overlap: max(min_A, min_B) <= min(max_A, max_B)
            lower = max(g_stats["ci_min"], b_stats["ci_min"])
            upper = min(g_stats["ci_max"], b_stats["ci_max"])
            
            is_overlap = lower <= upper
            res["Is_Overlap"] = is_overlap
            
            if is_overlap:
                to_keep_indices.add(i)
                if first_overlap_index == -1:
                    first_overlap_index = i
        
        # Transition Point Logic:
        # If we found an overlap at index N, we also want to keep N-1 (Result just before overlap).
        if first_overlap_index > 0:
             to_keep_indices.add(first_overlap_index - 1)
             print(f" -> Transition Point detected at Alpha {results[first_overlap_index - 1]['Alpha']}")
        elif first_overlap_index == 0:
             print(" -> Overlap starts immediately from highest alpha.")
        else:
             print(" -> No overlap detected in any dataset.")
             # Keep at least one result (the last one) if no overlap found
             if results:
                to_keep_indices.add(len(results) - 1)
             
        # Execute Deletion
        final_plot_data = []
        
        for i, res in enumerate(results):
            if i in to_keep_indices:
                final_plot_data.append(res)
                print(f" [Keep] Alpha {res['Alpha']} (Overlap: {res.get('Is_Overlap', 'N/A')})")
            else:
                print(f" [Delete] Alpha {res['Alpha']} (Stable, Redundant)")
                try:
                    shutil.rmtree(res["Folder"])
                except Exception as e:
                    print(f"    Error deleting {res['Folder']}: {e}")

        # Plotting
        self.plot_curve(final_plot_data, first_overlap_index=first_overlap_index, all_results=results)

    def plot_curve(self, data, first_overlap_index=-1, all_results=None):
        """
        7. Plots Rank G vs B with CIs.
        Generates the stability analysis plot and saves it to the results directory.
        """
        if not data:
            print("[Info] No data to plot.")
            return
            
        # Convert to DataFrame for easier plotting
        rows = []
        for d in data:
            rank_map = d["Stats"]
            row = {
                "Alpha": d["Alpha"] * 100, # Percent
                "Rank_G": rank_map.get("Method G", {}).get("rank"),
                "CI_Min_G": rank_map.get("Method G", {}).get("ci_min"),
                "CI_Max_G": rank_map.get("Method G", {}).get("ci_max"),
                
                "Rank_B": rank_map.get("Method B", {}).get("rank"),
                "CI_Min_B": rank_map.get("Method B", {}).get("ci_min"),
                "CI_Max_B": rank_map.get("Method B", {}).get("ci_max"),
                
                "Rank_D": rank_map.get("Method D", {}).get("rank"),
                "CI_Min_D": rank_map.get("Method D", {}).get("ci_min"),
                "CI_Max_D": rank_map.get("Method D", {}).get("ci_max"),
            }
            rows.append(row)
            
        df = pd.DataFrame(rows).sort_values(by="Alpha", ascending=True)
        
        # Plot Method G
        plt.plot(df["Alpha"], df["Rank_G"], label="Method G", color="blue", marker="o")
        plt.fill_between(df["Alpha"], df["CI_Min_G"], df["CI_Max_G"], color="blue", alpha=0.2)
        
        # Plot Method B
        plt.plot(df["Alpha"], df["Rank_B"], label="Method B", color="green", marker="x", linestyle="--")
        plt.fill_between(df["Alpha"], df["CI_Min_B"], df["CI_Max_B"], color="green", alpha=0.2)
        
        # Plot Method D
        if "Rank_D" in df.columns:
             plt.plot(df["Alpha"], df["Rank_D"], label="Method D", color="red", marker="^", linestyle=":")
             plt.fill_between(df["Alpha"], df["CI_Min_D"], df["CI_Max_D"], color="red", alpha=0.15)
             
        # Optimization: Add Overlap Line and Alpha Explanation
        if all_results and first_overlap_index >= 0:
             # Get the alpha value where overlap STARTS (the first index in the original sorted list that overlapped)
             # Note: results are sorted High Alpha -> Low Alpha. 
             # first_overlap_index is the index in that sorted list.
             # We want to mark the point where we transition from "No Overlap" to "Overlap".
             
             transition_alpha = all_results[first_overlap_index]["Alpha"] * 100
             
             plt.axvline(x=transition_alpha, color='black', linestyle='-.', alpha=0.7, label="Overlap Start")
             
             # Add Alpha Explanation Box
             textstr = '\n'.join((
                r'$\bf{Alpha}$: Contrast Level',
                r'(Higher values = Stronger Signal)',
                r' Vertical Line: Methods Indistinguishable'
                ))
             
             # Place a text box in upper right in axes coords
             props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
             plt.text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
                    verticalalignment='top', bbox=props)

        plt.xlabel("Alpha (%)")
        plt.ylabel("Rank (Lower is Better)")
        plt.title("Stability Analysis with Confidence Intervals")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.gca().invert_yaxis()
        
        save_path = self.results_dir / "stability_curve.png"
        plt.savefig(save_path, dpi=300)
        print(f" -> Optimization Curve saved to: {save_path}")

if __name__ == "__main__":
    workflow = HERAWorkflow()
    workflow.process_data()
