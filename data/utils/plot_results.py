#!/usr/bin/env python3
"""
Root CLI Entry Point for HERA Publication Plotting & Reporting Suite.

Syntax:
    python3 plot_results.py [target] [-o GRAPHICS_DIR] [-r REPORTS_DIR]

Description:
    Acts as the top-level convenience CLI wrapper forwarding execution to the
    encapsulated simulation plotting suite (utils/simulation/plot_results.py).
    Re-exports all plotting functions for 100% backwards compatibility with
    external scripts and publication workflows.

Author:
    Lukas von Erdmannsdorff
"""

import sys
from pathlib import Path

# Ensure the utils directory is on sys.path so 'simulation' is discoverable as a package
BASE_DIR = Path(__file__).parent.resolve()
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

from simulation.plot_results import *
from simulation.plot_results import __all__

if __name__ == "__main__":
    main()
