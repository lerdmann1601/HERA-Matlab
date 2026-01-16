"""
Prepares the HERA-Matlab Python distribution for PyPI publication.

This script patches the auto-generated MATLAB ``setup.py`` and ``__init__.py`` files
to ensure they meet PyPI metadata standards and include runtime safety checks.
It is intended to be run as part of the CI/CD pipeline.

Typical usage example:

    python .github/scripts/prepare_pypi.py [dist_dir]
"""

import glob
import os
import shutil
import sys
from typing import Optional


def prepare_distribution(target_dir: Optional[str] = None) -> None:
    """
    Main entry point for preparing the distribution.

    Locates the distribution directory, injects necessary assets (README, LICENSE),
    patches ``setup.py`` with enhanced metadata, and modifies ``__init__.py`` to
    include runtime availability checks.

    Args:
        target_dir: Optional path to the distribution directory. If not provided,
            the script attempts to auto-discover the directory.

    Raises:
        SystemExit: If the distribution directory cannot be found or if required
            files (like setup.py) are missing.
    """
    print("Pre-PyPI Distribution Preparation")

    # -------------------------------------------------------------------------
    # 1. Distribution Discovery
    # -------------------------------------------------------------------------
    dist_dir = ""
    
    if target_dir:
        if os.path.isdir(target_dir):
            dist_dir = target_dir
        else:
            print(f"Error: Provided directory '{target_dir}' does not exist.")
            sys.exit(1)
    else:
        # Attempt to find the distribution directory using glob patterns or default paths.
        dirs = [d for d in glob.glob("hera_matlab*") if os.path.isdir(d)]
        deploy_dir = "deploy/output/python"
        
        if dirs:
             dist_dir = dirs[0]
        elif os.path.isdir(deploy_dir) and os.path.exists(os.path.join(deploy_dir, "setup.py")):
             dist_dir = deploy_dir
        else:
            print("Error: No distribution directory found.")
            print(f"Current working directory: {os.getcwd()}")
            print(f"Contents: {os.listdir('.')}")
            sys.exit(1)
            
    print(f"Target Directory: {dist_dir}")

    # -------------------------------------------------------------------------
    # 2. Asset Integration
    # -------------------------------------------------------------------------
    print("Integrating Assets...")
    
    # Copy standard project files to the distribution root.
    files_to_copy = ["README.md", "LICENSE", "CITATION.cff"]
    for filename in files_to_copy:
        if os.path.exists(filename):
            print(f"  - Copying {filename}...")
            shutil.copy(filename, dist_dir)
        else:
            print(f"  - Warning: {filename} not found in root.")

    # Inject the runtime helper script into the package.
    runtime_script_source = "deploy/python_assets/install_runtime.py"
    if os.path.exists(runtime_script_source):
        # Locate the inner package directory (e.g., dist_dir/hera_matlab).
        pkg_inner_dir = os.path.join(dist_dir, "hera_matlab")
        if not os.path.isdir(pkg_inner_dir):
             print(f"  - Error: Package directory {pkg_inner_dir} not found.")
        else:
            runtime_script_dest = os.path.join(pkg_inner_dir, "install_runtime.py")
            print(f"  - Injecting {runtime_script_source} -> {runtime_script_dest}")
            shutil.copy(runtime_script_source, runtime_script_dest)
    else:
        print(f"  - Warning: {runtime_script_source} missing. Runtime check CLI will be unavailable.")

    # -------------------------------------------------------------------------
    # 3. Metadata Injection
    # -------------------------------------------------------------------------
    print("Patching setup.py...")
    setup_path = os.path.join(dist_dir, "setup.py")
    if not os.path.exists(setup_path):
        print("  - Error: setup.py not found.")
        sys.exit(1)

    with open(setup_path, "r") as f:
        content = f.read()

    # Code block to be injected into setup.py to enhance metadata.
    
    # 3a. Force Package Name Correction
    # Even if verified/patched before, we ensure the name is 'hera-matlab'
    # to match the PyPI Pending Publisher configuration.
    if "hera_matlab-R2025b" in content:
        print("  - Enforcing correct package name 'hera-matlab'...")
        content = content.replace("hera_matlab-R2025b", "hera-matlab")

    # 3b. Sync Version with GitHub Tag
    # MATLAB defaults to '25.2' (R2025b).
    # We overwrite this with the actual release tag (e.g., v1.1.0 -> 1.1.0).
    tag_name = os.environ.get('GITHUB_REF_NAME')
    if tag_name and tag_name.startswith('v'):
        new_version = tag_name[1:] # Strip 'v'
        print(f"  - Syncing version to GitHub Tag: {new_version}")
        import re
        # Replace version='25.2' with version='1.1.0'
        content = re.sub(r"version\s*=\s*['\"][\d\.]+['\"]", f"version='{new_version}'", content)

    with open(setup_path, "w") as f:
        f.write(content)
    injection_code = """
    # --- INJECTED METADATA START ---
    try:
        import os
        here = os.path.abspath(os.path.dirname(__file__))
        
        # Load Long Description from README.md
        readme_path = os.path.join(here, 'README.md')
        if os.path.exists(readme_path):
            with open(readme_path, 'r', encoding='utf-8') as f:
                setup_dict['long_description'] = f.read()
            setup_dict['long_description_content_type'] = 'text/markdown'
            
    except Exception as e:
        print(f'Warning: Could not inject README: {e}')

    # Enhanced Metadata
    setup_dict['name'] = 'hera-matlab'
    setup_dict['author'] = 'Lukas von Erdmannsdorff'
    setup_dict['license'] = 'MIT'
    
    # Keywords for PyPI discovery
    setup_dict['keywords'] = [
        'ranking', 'benchmarking', 'statistics', 'effect-size', 
        'bootstrapping', 'significance-testing', 'matlab-interface', 
        'hierarchical-compensatory', 'scientific-computing'
    ]

    # Dynamic URL handling based on GitHub Actions environment
    repo_url = 'https://github.com/lerdmann1601/HERA-Matlab' # Fallback
    if 'GITHUB_REPOSITORY' in os.environ:
         repo_url = f"https://github.com/{os.environ['GITHUB_REPOSITORY']}"
    
    setup_dict['url'] = repo_url
    setup_dict['project_urls'] = {
        'Bug Tracker': f"{repo_url}/issues",
        'Source Code': repo_url,
        'Documentation': f"{repo_url}/tree/main/docs",
    }
    
    setup_dict['classifiers'] = [
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]

    # Python Version Constraint
    setup_dict['python_requires'] = '>=3.9, <3.13'
    # --- INJECTED METADATA END ---
    
    setup(**setup_dict)
    """

    if "INJECTED METADATA START" in content:
        print("  - Metadata already injected. Skipping setup.py patch.")
    elif "setup(**setup_dict)" in content:
        new_content = content.replace("setup(**setup_dict)", injection_code)
        with open(setup_path, "w") as f:
            f.write(new_content)
        print("  - Success: setup.py patched.")
    else:
        print("  - Error: Could not find anchor point 'setup(**setup_dict)' in setup.py.")

    # -------------------------------------------------------------------------
    # 4. Runtime Safety Injection
    # -------------------------------------------------------------------------
    print("Patching __init__.py for Runtime Safety...")
    init_path = os.path.join(dist_dir, "hera_matlab", "__init__.py")
    if os.path.exists(init_path):
        with open(init_path, "r") as f:
            init_content = f.read()
            
        if "INJECTED RUNTIME CHECK" in init_content:
             print("  - Runtime check already injected. Skipping __init__.py patch.")
        else:
            # The pattern generated by MATLAB's setup.
            strict_init_pattern = """_pir = _PathInitializer()
_pir.get_paths_from_os()
_pir.update_paths()
_pir.import_cppext()
_pir.import_matlab_pysdk_runtime()
_pir.import_matlab()"""

            # The robust block that catches initialization errors.
            robust_init_block = """
# --- INJECTED RUNTIME CHECK ---
_initialization_error = None
try:
    _pir = _PathInitializer()
    _pir.get_paths_from_os()
    _pir.update_paths()
    _pir.import_cppext()
    _pir.import_matlab_pysdk_runtime()
    _pir.import_matlab()
except (RuntimeError, ImportError, EnvironmentError) as e:
    _initialization_error = e
    pass
# --- INJECTED RUNTIME CHECK END ---
"""
            
            init_func_pattern = "def initialize():"
            init_check_code = """def initialize():
    if _initialization_error:
        print(f"Error: Failed to initialize MATLAB Runtime.")
        print(f"Details: {str(_initialization_error)}")
        print("")
        print("Run the following command to diagnose and fix the issue:")
        print("    python -m hera_matlab.install_runtime")
        raise _initialization_error
"""

            if strict_init_pattern in init_content:
                init_content = init_content.replace(strict_init_pattern, robust_init_block)
                init_content = init_content.replace(init_func_pattern, init_check_code)
                with open(init_path, "w") as f:
                    f.write(init_content)
                print("  - Success: __init__.py patched.")
            else:
                 print("  - Warning: Could not find strict initialization block in __init__.py.")
    else:
        print(f"  - Error: {init_path} not found.")

    print("Preparation Complete.")

if __name__ == "__main__":
    target = None
    if len(sys.argv) > 1:
        target = sys.argv[1]
    prepare_distribution(target)
