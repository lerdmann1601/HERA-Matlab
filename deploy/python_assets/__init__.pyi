"""
Type stubs for HERA-Matlab Python Package.
This file enables IntelliSense/Autocompletion for Python users.
"""

from typing import Optional, Dict, Any, Union

def start_ranking(configFile: Optional[str] = ..., runtest: Optional[str] = ..., logPath: Optional[str] = ...) -> None:
    """
    Interactive script to configure and start the ranking.
    
    Args:
        configFile: Path to a JSON configuration file for batch mode.
        runtest: Set to 'true' to run unit tests.
        logPath: Path to save unit test logs.
    """
    ...

def run_ranking(userInput: Dict[str, Any]) -> None:
    """
    Main entry point for running the HERA ranking analysis.
    
    Args:
        userInput: Configuration dictionary containing the analysis parameters.
        
    Example:
    Example:
        # Direct Injection Mode (NumPy arrays are automatically converted)
        config = {
            "userInput": {
                "custom_data": [numpy_array_1, numpy_array_2],
                "metric_names": ["Runtime", "Accuracy"],
                "ranking_mode": "M1_M2",
                # ... other options
            }
        }
        hera.run_ranking(config)
    """
    ...

def initialize() -> None:
    """
    Initialize the MATLAB Runtime.
    Raises an error if the runtime cannot be found.
    """
    ...

def terminate() -> None:
    """
    Terminate the MATLAB Runtime.
    """
    ...
