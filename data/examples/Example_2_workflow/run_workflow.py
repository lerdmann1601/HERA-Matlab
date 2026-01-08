"""
RUN_WORKFLOW - Example workflow for HERA analysis with pruning and stability checking.

Description:
  This script demonstrates a complete workflow for running HERA on a series of datasets 
  (e.g., decreasing contrast/alpha levels). It performs the following steps:
  1. Sets up the necessary directory structure.
  2. Iterates through data folders (sorted by alpha).
  3. Runs the HERA analysis (via compiled runtime or local script).
  4. Analyzes the results, calculating Confidence Intervals (CIs) from bootstrap ranks of not using pre-calculated CIs from JSON.
  5. Prunes the results to retain only the transition point and stable regions where methods overlap.
  6. Generates a stability curve plot.

Workflow:
  1. Initialization: Define paths and configuration.
  2. Setup: Ensure data and results directories exist.
  3. Data Processing: 
     - Detect input folders.
     - Configure HERA for each dataset.
     - Execute HERA.
     - Parse JSON output.
  4. Pruning & Plotting:
     - Identify overlap regions between methods (e.g., Method G and Method B).
     - Keep relevant datasets (overlapping + transition point).
     - Remove redundant stable datasets.
     - Plot the optimization curve.

Usage:
  Run this script directly from the terminal or IDE:
  $ pip install -r requirements.txt
  $ python3 run_workflow.py

Author: Lukas von Erdmannsdorff
"""

import os
import time
import sys
import shutil
import json
import subprocess
import platform
import glob
from pathlib import Path
from typing import Any

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


# --- Configuration ---
BASE_DIR = Path(__file__).parent.resolve()
DATA_DIR = BASE_DIR / "data"
RESULTS_DIR = BASE_DIR / "results"
# Determine OS-specific defaults
current_os = platform.system()

if current_os == "Darwin":
    DEFAULT_DEPLOY_DIR = BASE_DIR.parent.parent.parent / "release" / "macos"
    DEFAULT_EXE_NAME = "run_HERA_Runtime.sh"
    # Default MCR path for Mac (User should verify this!)
    DEFAULT_MCR_PATH = "/Applications/MATLAB_R2025b.app" 
elif current_os == "Windows":
    DEFAULT_DEPLOY_DIR = BASE_DIR.parent.parent.parent / "release" / "windows"
    DEFAULT_EXE_NAME = "HERA_Runtime.exe" # Verify this name
    DEFAULT_MCR_PATH = "C:\\Program Files\\MATLAB\\R2025b"
elif current_os == "Linux":
    DEFAULT_DEPLOY_DIR = BASE_DIR.parent.parent.parent / "release" / "linux"
    DEFAULT_EXE_NAME = "run_HERA_Runtime.sh"
    DEFAULT_MCR_PATH = "/usr/local/MATLAB/R2025b"
else:
    DEFAULT_DEPLOY_DIR = BASE_DIR.parent.parent.parent / "release"
    DEFAULT_EXE_NAME = "run_HERA_Runtime.sh"
    DEFAULT_MCR_PATH = ""

HERA_DEPLOY_DIR_ENV = os.environ.get("HERA_DEPLOY_DIR")
HERA_EXECUTABLE_NAME_ENV = os.environ.get("HERA_EXECUTABLE_NAME")
HERA_MCR_PATH_ENV = os.environ.get("HERA_MCR_PATH")

METRIC_NAMES = ["OC", "CNR", "SNR"]
RANKING_MODE = "M1_M2_M3"

# Manual Bootstrap Settings
MANUAL_BOOTSTRAP = {
    "manual_B_thr": 2000,
    "manual_B_ci": 5000,
    "manual_B_rank": 500
}

class Logger:
    def __init__(self, file_path: Path | str) -> None:
        self.terminal = sys.stdout
        self.log = open(file_path, "a")

    def write(self, message: str) -> None:
        self.terminal.write(message)
        self.log.write(message)
        self.log.flush()  

    def flush(self):
        # needed for python 3 compatibility
        self.terminal.flush()
        self.log.flush()

class HERAWorkflow:
    """
    Main controller for the HERA analysis workflow.
    
    This class manages the lifecycle of the analysis, including directory setup, configuration generation,
    execution of the HERA runtime, and result processing.
    """

    def __init__(self) -> None:
        """Initializes the workflow with default data and results directories."""
        self.data_dir: Path = DATA_DIR
        self.results_dir: Path = RESULTS_DIR
        
    def setup_directories(self) -> None:
        """
        1. Ensures the directory structure exists.
        
        Checks for data and result directories and creates them if missing.
        """
        # 1. Setup Data Directory
        if not self.data_dir.exists():
            print(f"[Info] Creating Data Directory: {self.data_dir}")
            self.data_dir.mkdir(parents=True)
            
        # 2. Setup Results Directory
        if not self.results_dir.exists():
            print(f"[Info] Creating Results Directory: {self.results_dir}")
            self.results_dir.mkdir(parents=True)

        # --- Logging Setup ---
        # Generate timestamp
        timestamp = time.strftime("%Y%m%d_%H%M%S")
        log_file_name = f"workflow_log_{timestamp}.txt"
        log_path = self.results_dir / log_file_name
        
        # Redirect stdout
        sys.stdout = Logger(log_path)
        print(f"[Info] Logging started. Output saved to: {log_path}")
                 
    def create_hera_config(self, input_path: Path, output_path: Path) -> dict[str, Any]:
        """
        2. Creates the JSON configuration structure for HERA.
        
        Constructs the dictionary expected by the HERA runtime.

        Args:
            input_path (Path): Path to the input data directory.
            output_path (Path): Path where results should be saved.

        Returns:
            dict[str, Any]: A dictionary containing the 'userInput' configuration.
        """
        config = {
            "userInput": {
                # --- Mandatory ---
                "folderPath": str(input_path),
                "metric_names": METRIC_NAMES,
                "output_dir": str(output_path), # Output directory for results would be "results/Alpha_XX" otherwise by default..
                
                # --- Optional ---
                "fileType": ".csv",
                "create_reports": False,  # Disabled for speed
                
                # --- Manual Bootstrap ---
                "manual_B_thr": MANUAL_BOOTSTRAP["manual_B_thr"],
                "manual_B_ci": MANUAL_BOOTSTRAP["manual_B_ci"],
                "manual_B_rank": MANUAL_BOOTSTRAP["manual_B_rank"]
            }
        }
        return config

    def find_executable(self) -> Path | str | None:
        """
        Helper: Searches for the HERA executable in likely locations.
        
        Returns:
            Path | str | None: The path to the executable (Path object), "matlab" (str) if using local MATLAB, 
                                        or None if not found.
        """
        # 1. Check Environment Variable
        if HERA_DEPLOY_DIR_ENV:
             custom_path = Path(HERA_DEPLOY_DIR_ENV) / (HERA_EXECUTABLE_NAME_ENV or DEFAULT_EXE_NAME)
             if custom_path.exists():
                 return custom_path
        
        # 2. Check Release/Deploy and Standard Installation Paths
        potential_dirs = [
            BASE_DIR.parent.parent.parent / "release" / "macos",   # Release structure (Mac)
            BASE_DIR.parent.parent.parent / "release" / "windows", # Release structure (Win)
            BASE_DIR.parent.parent.parent / "release" / "linux",   # Release structure (Linux)
            BASE_DIR.parent.parent.parent / "deploy" / "output",   # Local Build
            Path("/Applications/HERA_Runtime/application"),               # MacOS Default
            Path("C:/Program Files/HERA_Runtime/application"),            # Windows Default
            Path(f"{Path.home()}/Applications/HERA_Runtime/application"), # MacOS/Linux User Local
            Path("/usr/local/HERA_Runtime/application"),                  # Linux System Default
        ]
        
        target_name = HERA_EXECUTABLE_NAME_ENV or DEFAULT_EXE_NAME
        
        for p_dir in potential_dirs:
            exe_path = p_dir / target_name
            if exe_path.exists():
                return exe_path
                
        # 3. Fallback: HERA_Launcher.command (Mac specific convenience)
        launcher_path = BASE_DIR.parent.parent.parent / "release" / "macos" / "HERA_Launcher.command"
        if launcher_path.exists():
            return launcher_path

        # 4. Fallback: Check for MATLAB executable in PATH
        matlab_path = shutil.which("matlab")
        if matlab_path:
             return "matlab"

        # 5. Fallback: Check for local MATLAB installation (MacOS)
        if current_os == "Darwin":
             # Look for standard MATLAB installations
             candidates = sorted(glob.glob("/Applications/MATLAB_*.app/bin/matlab"))
             if candidates:
                 # Use the latest version (last alphabetically usually works for R20xx)
                 return candidates[-1]

        return None

    def run_hera(self, config_path: Path) -> bool:
        """
        3. Executes HERA using the compiled runtime/launcher.
        
        Locates the executable and runs it with the provided configuration.

        Args:
            config_path (Path): Path to the generated JSON configuration file.
            
        Returns:
            bool: True if execution was successful, False otherwise.
        """
        hera_cmd = self.find_executable()
        
        mcr_path = HERA_MCR_PATH_ENV or DEFAULT_MCR_PATH

        if not hera_cmd:
             print(f"\n[Error] HERA executable NOT found and MATLAB is NOT in the system path!")
             print(f"  To fix this, either:")
             print(f"  1. Build the project using 'deploy/build_hera.m'")
             print(f"  2. Download the release from GitHub")
             print(f"  3. Ensure 'matlab' is accessible from the terminal (add to PATH)")
             print(f"  Alternatively, set the 'HERA_DEPLOY_DIR' environment variable.\n")
             return False

        print(f" -> Running HERA using: {hera_cmd}")
        print(f" -> Config: {config_path}")
        
        # Check if we are using MATLAB (System path or local path)
        is_matlab = str(hera_cmd).endswith("matlab") or str(hera_cmd) == "matlab"

        # Branch 1: Run via Compiled Executable
        if not is_matlab:
            cmd = [str(hera_cmd), mcr_path, "configFile", str(config_path)]
            try:
                result = subprocess.run(cmd, capture_output=True, text=True)
                if result.returncode != 0:
                    print("[Error] HERA Runtime Failed!")
                    print(result.stderr)
                    # print(result.stdout) # optional debug
                    return False
                return True
            except Exception as e:
                print(f"[Error] Execution Exception: {e}")
                return False
                
        # Branch 2: Run via MATLAB Command Line (No compiled binary needed)
        else:
            # We must construct a MATLAB batch command that:
            # 1. Adds HERA to path (using setup_HERA in the repo root)
            # 2. Runs start_ranking with the config
            
            repo_root = BASE_DIR.parent.parent.parent.resolve()
            
            # Note: We need to cd to repo_root to ensure setup_HERA is found, 
            # OR add repo_root to path explicitly in the command.
            # Simpler: cd to repo_root in the subprocess call.
            
            matlab_cmd_str = f"setup_HERA; HERA.start_ranking('configFile', '{config_path}');"
            
            # Use the found hera_cmd (which might be a full path to matlab binary)
            cmd = [str(hera_cmd), "-batch", matlab_cmd_str]
            
            print(f" -> Executing MATLAB: {matlab_cmd_str}")
            
            try:
                # Run in repo_root so "setup_HERA" is found immediately
                result = subprocess.run(cmd, cwd=repo_root, capture_output=True, text=True)
                if result.returncode != 0:
                     print("[Error] MATLAB Execution Failed!")
                     print(result.stderr)
                     print(result.stdout)
                     return False
                return True
            except Exception as e:
                print(f"[Error] MATLAB subprocess error: {e}")
                return False

    def get_latest_result_json(self, result_folder: Path) -> Path | None:
        """
        Helper: Finds the latest JSON output in the results folder.

        Args:
            result_folder (Path): The folder containing HERA results.

        Returns:
            Path | None: Path to the latest result JSON file, or None if not found.
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

    def analyze_results(self, json_path: Path) -> dict[str, Any]:
        """
        4. Analysis: Reads ranking info and calculates CIs from bootstrap ranks.
        
        Specific logic to extract bootstrap distributions and compute percentiles.

        Args:
            json_path (Path): Path to the HERA result JSON.

        Returns:
            dict[str, Any]: A dictionary containing rank and CI statistics for each method.
        """
        with open(json_path, 'r') as f:
            data = json.load(f)
            
        results: dict[str, Any] = data.get('results', {})
        dataset_names: list[str] = data.get('dataset_names', [])
        
        final_ranks = results.get('final_rank', [])
        
        # Calculate CIs from bootstrap distribution or use pre-calculated CIs
        ci_lower = results.get('ci_lower_rank', [])
        ci_upper = results.get('ci_upper_rank', [])
        
        if ci_lower and ci_upper and len(ci_lower) == len(dataset_names):
             print(f" -> Using pre-calculated CIs from JSON.")
        else:
            # Fallback: Calculate from bootstrap ranks
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

        rank_map: dict[str, dict[str, Any]] = {}
        for i, name in enumerate(dataset_names):
            rank_map[name] = {
                "rank": final_ranks[i],
                "ci_min": ci_lower[i],
                "ci_max": ci_upper[i]
            }
            
        return rank_map

    def process_data(self) -> None:
        """
        Main Loop: Iterates through data folders, runs HERA, and collects results.
        
        This method orchestrates the entire batch processing workflow:
        1. directory setup
        2. processing loop (config -> run -> analyze)
        3. result aggregation
        4. handing off to pruning/plotting.
        """
        self.setup_directories()
 
        # Find all Alpha_XX folders in python_experiments/data
        alpha_folders = [d for d in self.data_dir.iterdir() if d.is_dir() and d.name.startswith("Alpha_")]
        
        # Sort by alpha value descending (High Alpha/Contrast -> Low Alpha/Noise)
        def parse_alpha(p):
            try:
                return int(p.name.split("_")[-1])
            except (ValueError, IndexError):
                return -1
                
        alpha_folders.sort(key=parse_alpha, reverse=True)
        
        if not alpha_folders:
            print("[Info] No data folders found to process.")
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

    def prune_and_plot(self, results: list[dict[str, Any]]) -> None:
        """
        5. Pruning and Plotting Logic.
        
        Applies pruning rules to the aggregated results and triggers plotting.
        
        Pruning Rules:
          1. Keep all results where Method G and Method B CIs OVERLAP.
          2. Keep the method immediately BEFORE the overlap starts (Transition point).
          3. Delete all others (Stable results).
          
        Args:
            results (list[dict[str, Any]]): List of result dictionaries from `process_data`.
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

    def plot_curve(self, data: list[dict[str, Any]], first_overlap_index: int = -1, all_results: list[dict[str, Any]] | None = None) -> None:
        """
        6. Plots Rank G vs B with CIs.
        
        Generates the stability analysis plot and saves it to the results directory.
        
        Args:
            data (list[dict[str, Any]]): Filtered list of results to be plotted.
            first_overlap_index (int): Index where overlap was first detected (for vertical line).
            all_results (list[dict[str, Any]] | None): Full list of results for transition point lookup.
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
             
        # Add Overlap Line and Alpha Explanation
        if all_results and first_overlap_index >= 0:
             # Get the alpha value where overlap STARTS (the first index in the original sorted list that overlapped)
             # Note: results are sorted High Alpha -> Low Alpha. 
             # first_overlap_index is the index in that sorted list.
             # We want to mark the point where we transition from "No Overlap" to "Overlap".
             
             transition_alpha = all_results[first_overlap_index]["Alpha"] * 100
             
             plt.axvline(x=transition_alpha, color='black', linestyle='-.', alpha=0.7, label="CI-Overlap")
             
             # Add Alpha Explanation Box (Below the plot)
             textstr = '\n'.join((
                r'$\bf{Alpha}$: Enhancement Level (Higher values = Stronger Contrast)',
                ))
             
             # Adjust subplot to make room at bottom
             plt.subplots_adjust(bottom=0.25)
             
             # Place text below x-axis
             props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
             plt.figtext(0.5, 0.02, textstr, ha='center', fontsize=10, bbox=props)

        plt.xlabel("Alpha (%)")
        plt.ylabel("Rank (Lower is Better)")
        plt.title("Stability Analysis using Confidence Intervals")
        plt.legend()
        plt.grid(True, alpha=0.3)
        plt.gca().invert_yaxis()
        
        save_path = self.results_dir / "stability_curve.png"
        plt.savefig(save_path, dpi=300)
        print(f" -> Optimization Curve saved to: {save_path}")

if __name__ == "__main__":
    start_time = time.time()
    workflow = HERAWorkflow()
    workflow.process_data()
    end_time = time.time()
    duration = end_time - start_time
    
    if duration < 60:
        print(f"Total analysis duration: {duration:.2f} seconds")
    elif duration < 3600:
        minutes = int(duration // 60)
        seconds = duration % 60
        print(f"Total analysis duration: {minutes} min {seconds:.0f} sec")
    else:
        hours = int(duration // 3600)
        minutes = int((duration % 3600) // 60)
        seconds = duration % 60
        print(f"Total analysis duration: {hours} h {minutes} min {seconds:.0f} sec")
