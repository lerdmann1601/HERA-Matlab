"""
PREPARE_PYPI - Prepares the Python distribution for PyPI publication.

Syntax:
  python .github/scripts/prepare_pypi.py

Description:
  This script is part of the Continuous Integration (CI) pipeline for publishing to PyPI.
  It works around limitations of the auto-generated MATLAB setup.py by:
  1.  Injecting the correct project 'README.md' as the long description.
  2.  Setting the correct 'url' metadata to point to the GitHub repository.
  3.  Ensuring all metadata is compliant with modern PyPI standards.

Workflow:
  1.  Distribution Discovery: Find the extracted source distribution directory.
  2.  Asset Integration: Copy 'README.md' from the repository root into the package source.
  3.  Metadata Injection: Patch 'setup.py' to dynamicall load the README and repository URL.

Author: Lukas von Erdmannsdorff
"""

import os
import sys
import shutil
import glob

def prepare_distribution():
    # 1. Distribution Discovery
    # Find the directory extracted from the source tarball (e.g., 'hera_matlab-25.2/').
    dirs = [d for d in glob.glob("hera_matlab-*/") if os.path.isdir(d)]
    if not dirs:
        print("Error: No distribution directory found.")
        sys.exit(1)
    
    dist_dir = dirs[0]
    print(f"Found distribution directory: {dist_dir}")

    # 2. Asset Integration
    # Copy main project files (README.md, LICENSE) into the distribution folder.
    files_to_copy = ["README.md", "LICENSE"]
    for filename in files_to_copy:
        if os.path.exists(filename):
            print(f"Copying {filename}...")
            shutil.copy(filename, dist_dir)
        else:
            print(f"Warning: No {filename} found in root.")

    # 3. Metadata Injection
    # Locate the setup.py file that needs to be modified.
    setup_path = os.path.join(dist_dir, "setup.py")
    if not os.path.exists(setup_path):
        print("Error: setup.py not found.")
        sys.exit(1)

    print(f"Patching {setup_path}...")
    with open(setup_path, "r") as f:
        content = f.read()

    # Define the Python code block to be injected into setup.py.
    # This code runs at install time to load the README and set the URL.
    injection_code = """
    # --- INJECTED METADATA START ---
    try:
        import os
        here = os.path.abspath(os.path.dirname(__file__))
        
        # Load Long Description
        readme_path = os.path.join(here, 'README.md')
        if os.path.exists(readme_path):
            with open(readme_path, 'r', encoding='utf-8') as f:
                setup_dict['long_description'] = f.read()
            setup_dict['long_description_content_type'] = 'text/markdown'
            
    except Exception as e:
        print(f'Warning: Could not inject README: {e}')

    # Enhanced Metadata
    setup_dict['author'] = 'Lukas von Erdmannsdorff'
    setup_dict['license'] = 'MIT'
    
    # Keywords for discovery
    setup_dict['keywords'] = [
        'ranking', 'benchmarking', 'statistics', 'effect-size', 
        'bootstrapping', 'significance-testing', 'matlab-interface', 
        'hierarchical-compensatory', 'scientific-computing'
    ]

    # Dynamic URL handling based on CI environment
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
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Mathematics',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Libraries :: Python Modules',
    ]
    # --- INJECTED METADATA END ---
    
    setup(**setup_dict)
    """

    # Apply the patch by replacing the original setup() call with our injected logic.
    if "setup(**setup_dict)" in content:
        new_content = content.replace("setup(**setup_dict)", injection_code)
        
        with open(setup_path, "w") as f:
            f.write(new_content)
        print("Successfully patched setup.py")
    else:
        print("Error: Could not find 'setup(**setup_dict)' in setup.py to patch.")
        sys.exit(1)

if __name__ == "__main__":
    prepare_distribution()
