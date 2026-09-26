#!/usr/bin/env python3
"""
Root CLI Entry Point for HERA Methodological Ground-Truth Validation.

Syntax:
    python3 run_simulation.py [options]

Description:
    Acts as the top-level convenience CLI wrapper forwarding execution to the
    encapsulated simulation package (utils/simulation/). Re-exports all public
    symbols for 100% backwards compatibility with external scripts and workflows.

Author:
    Lukas von Erdmannsdorff
"""

import sys
from pathlib import Path

# Ensure the utils directory is on sys.path so 'simulation' is discoverable as a package
BASE_DIR = Path(__file__).parent.resolve()
if str(BASE_DIR) not in sys.path:
    sys.path.insert(0, str(BASE_DIR))

# Import and re-export all symbols from the encapsulated package for 100% backwards compatibility
from simulation.run_simulation import *
from simulation.run_simulation import __all__

if __name__ == "__main__":
    main()
