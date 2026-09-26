"""
HERA Computational Engine Wrapper.

Manages execution of the HERA ranking core across both the compiled Python package (PyPI)
and the local MATLAB CLI backend with isolated single-thread allocation (-singleCompThread).

Author: Lukas von Erdmannsdorff
"""

import sys
import glob
import shutil
import platform
import logging
import subprocess
from pathlib import Path
from typing import Optional

from .system import suppress_stdout_stderr
from .config import BASE_DIR


class HERAExecutor:
    """
    Manages execution of HERA across Python Package (PyPI) and local MATLAB environments.

    Syntax:
        executor = HERAExecutor(logger, repo_root)

    Description:
        Wraps HERA evaluation core. Detects whether the compiled binary package
        (`hera_matlab`) is present or falls back to invoking the local MATLAB CLI
        in headless `-nodisplay -batch -singleCompThread` mode.

    Parameters:
        logger (Optional[logging.Logger]): Logger instance for logging environment discovery.
        repo_root (Optional[Path]): Root path of the HERA-Matlab repository.

    Author:
        Lukas von Erdmannsdorff
    """

    def __init__(self, logger: Optional[logging.Logger] = None, repo_root: Optional[Path] = None):
        self.logger = logger
        self.engine = None
        self.mode = None  # "python_pkg" or "matlab_cli"
        self.matlab_bin = None
        self.repo_root = repo_root or BASE_DIR.parent.parent.parent.resolve()

        # 1. Attempt PyPI package initialization with suppressed native stdout
        with suppress_stdout_stderr():
            try:
                import hera_matlab
                try:
                    self.engine = hera_matlab.initialize()
                except AttributeError:
                    from hera_matlab.runtime_wrapper import HeraSmartWrapper
                    self.engine = HeraSmartWrapper(hera_matlab._pir.initialize_package())
                self.mode = "python_pkg"
            except Exception:
                self.engine = None

        if self.mode == "python_pkg":
            if self.logger:
                self.logger.info("[Environment] HERA computational engine initialized via Python package (hera-matlab).")
            return

        # 2. Fallback to MATLAB CLI
        matlab_cmd = shutil.which("matlab")
        if not matlab_cmd and platform.system() == "Darwin":
            candidates = sorted(glob.glob("/Applications/MATLAB_*.app/bin/matlab"))
            if candidates:
                matlab_cmd = candidates[-1]

        if matlab_cmd:
            self.matlab_bin = matlab_cmd
            self.mode = "matlab_cli"
            if self.logger:
                self.logger.info(f"[Environment] Found local MATLAB executable: {self.matlab_bin}")
        else:
            raise RuntimeError("Neither hera_matlab package nor local MATLAB executable could be found.")

    def run(self, config_path: Path, timeout: float = 300.0) -> None:
        """
        Executes HERA ranking using the active backend with isolated single-thread allocation.

        Syntax:
            executor.run(config_path, timeout=300.0)

        Parameters:
            config_path (Path): Path to the temporary JSON configuration file for the run.
            timeout (float): Execution timeout in seconds (default: 300.0).

        Returns:
            None

        Author:
            Lukas von Erdmannsdorff
        """
        with suppress_stdout_stderr():
            if self.mode == "python_pkg":
                self.engine.start_ranking("configFile", str(config_path.resolve()), nargout=0)
            elif self.mode == "matlab_cli":
                cmd_str = f"setup_HERA; HERA.start_ranking('configFile', '{config_path.resolve()}');"
                # Enforce -singleCompThread to prevent core oversubscription among parallel workers
                cmd = [
                    str(self.matlab_bin),
                    "-nodisplay",
                    "-nosplash",
                    "-singleCompThread",
                    "-batch",
                    cmd_str
                ]
                res = subprocess.run(cmd, cwd=self.repo_root, capture_output=True, text=True, timeout=timeout)
                if res.returncode != 0:
                    raise RuntimeError(f"MATLAB execution failed (code {res.returncode}):\n{res.stderr}\n{res.stdout}")

    def terminate(self) -> None:
        """
        Cleans up engine resources if running under Python package.

        Syntax:
            executor.terminate()

        Returns:
            None

        Author:
            Lukas von Erdmannsdorff
        """
        if self.mode == "python_pkg" and self.engine is not None:
            try:
                with suppress_stdout_stderr():
                    self.engine.terminate()
            except Exception:
                pass


# Global worker executor singleton for multiprocessing worker processes
_WORKER_EXECUTOR: Optional[HERAExecutor] = None


def get_worker_executor(repo_root: Path) -> HERAExecutor:
    """
    Returns or initializes worker-local HERA executor instance singleton.

    Syntax:
        executor = get_worker_executor(repo_root)

    Parameters:
        repo_root (Path): Root directory of the HERA-Matlab repository.

    Returns:
        HERAExecutor: Initialized executor instance for the current worker process.

    Author:
        Lukas von Erdmannsdorff
    """
    global _WORKER_EXECUTOR
    if _WORKER_EXECUTOR is None:
        _WORKER_EXECUTOR = HERAExecutor(logger=None, repo_root=repo_root)
    return _WORKER_EXECUTOR
