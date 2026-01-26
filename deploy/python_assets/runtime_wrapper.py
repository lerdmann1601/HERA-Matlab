
import sys
from typing import Any, Dict, List, Union

try:
    import numpy as np
    HAS_NUMPY = True
except ImportError:
    HAS_NUMPY = False

try:
    import matlab
    HAS_MATLAB = True
except ImportError:
    HAS_MATLAB = False

class HeraSmartWrapper:
    """
    A smart wrapper around the MATLAB Runtime instance.
    
    This class acts as a proxy, intercepting specific method calls to provide
    enhanced functionality (like automatic data conversion) while delegating
    all other calls to the underlying MATLAB Runtime object.

    Attributes:
        _hera: The original MATLAB Runtime instance.
    """
    
    def __init__(self, original_instance: Any):
        """
        Initializes the wrapper.

        Args:
            original_instance: The instantiated MATLAB Runtime object.
        """
        self._hera = original_instance

    def __getattr__(self, name: str) -> Any:
        """
        Delegates attribute access to the underlying instance.

        Args:
            name: The name of the attribute to access.

        Returns:
            The attribute from the original MATLAB Runtime instance.
        """
        return getattr(self._hera, name)

    def run_ranking(self, userInput: Dict[str, Any], **kwargs) -> Any:
        """
        Intercepts the run_ranking call to perform automatic type conversion.

        Recursively converts any NumPy arrays found within the input dictionary
        into MATLAB-compatible types (e.g., matlab.double) before passing them
        to the runtime.

        Args:
            userInput: The configuration dictionary for the ranking.
            **kwargs: Additional keyword arguments passed to the runtime.

        Returns:
            The result from the MATLAB run_ranking function.
        """
        if HAS_NUMPY and HAS_MATLAB:
            # We only attempt conversion if both libs are present
            userInput = self._convert_numpy_to_matlab(userInput)
        
        return self._hera.run_ranking(userInput, **kwargs)

    def _convert_numpy_to_matlab(self, data: Any) -> Any:
        """
        Recursively searches for NumPy arrays and converts them to matlab.double.
        """
        if isinstance(data, dict):
            return {k: self._convert_numpy_to_matlab(v) for k, v in data.items()}
        elif isinstance(data, list):
            return [self._convert_numpy_to_matlab(v) for v in data]
        elif isinstance(data, np.ndarray):
            # Check for numeric type
            if np.issubdtype(data.dtype, np.number):
                try:
                    return matlab.double(data.tolist())
                except Exception:
                    # Fallback or pass through if conversion fails
                    return data
            else:
                return data
        else:
            return data
