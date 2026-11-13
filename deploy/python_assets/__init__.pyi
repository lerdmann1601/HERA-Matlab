"""
Type stubs for HERA-Matlab Python Package.
This file enables IntelliSense/Autocompletion for Python users.
"""

from typing import Optional, Dict, Any, Union

__version__: str = ...

def start_ranking(
    configFile: Optional[str] = ...,
    runtest: Optional[str] = ...,
    logPath: Optional[str] = ...,
    convergence: Optional[str] = ...,
    sims: Optional[int] = ...,
    **kwargs: Any
) -> None:
    """
    Interactive script to configure and start the ranking, unit tests, or convergence analysis.
    
    Note: A version check is performed upon initialization to notify the user of updates.
    
    Args:
        configFile: Path to a JSON configuration file for batch mode.
        runtest: Set to 'true' to run unit tests.
        logPath: Path to save unit test or convergence logs.
        convergence: Set to 'true' or a JSON path to run the convergence study.
        sims: Number of simulations per scenario for convergence analysis.
        **kwargs: Additional MATLAB runtime arguments (e.g., nargout=0).
    """
    ...

def run_ranking(userInput: Dict[str, Any], **kwargs: Any) -> Dict[str, Any]:
    """
    Main entry point for running the HERA ranking analysis.
    
    Args:
        userInput: Configuration dictionary containing the analysis parameters.
        **kwargs: Additional MATLAB runtime arguments (e.g., nargout=1).
        
    Returns:
        Dictionary containing all computed statistics and ranking results.
        
    Example:
        # Direct Injection Mode (NumPy arrays are automatically converted)
        config = {
            "custom_data": [numpy_array_1, numpy_array_2],
            "metric_names": ["Throughput", "Accuracy"],
            "ranking_mode": "M1_M2",
            # ... other options
        }
        results = hera.run_ranking(config, nargout=1)
    """
    ...

def initialize() -> Any:
    """
    Initialize the MATLAB Runtime and return the HERA smart wrapper instance.
    Raises an error if the runtime cannot be found.
    """
    ...

def terminate() -> None:
    """
    Terminate the MATLAB Runtime.
    """
    ...
