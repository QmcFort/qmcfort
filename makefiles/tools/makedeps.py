import os
import re
import sys

# Analyze script's argument list and extract source files directories
# and corresponding build directories
if (len(sys.argv) == 3):
    src_dirs = [sys.argv[1]]
    build_dirs = [sys.argv[2]]
elif (len(sys.argv) == 5):
    src_dirs = [sys.argv[1], sys.argv[3]]
    build_dirs = [sys.argv[2], sys.argv[4]]
else:
    src_dirs = [os.getcwd()]
    build_dirs = [os.getcwd()]

# Output dependency file
DEPEND_FILE = "deps.mk"

# Regex pattern to match 'use' statements including cases with 'only'
use_pattern = re.compile(r"^\s*use\s+(\w+)(?:,\s*only\s*:\s*[\w\s,]*)?", re.IGNORECASE)

# Regex pattern to match 'module' statements
module_pattern = re.compile(r"^\s*module\s+(\w+)", re.IGNORECASE)

# List to store all Fortran files
file_list = []

# Store object and source paths for each of the files
obj_list = {}
src_list = {}

# Dictionary to store all user-defined modules and the files where they are defined
# module_dict[module_name] = fname
module_dict = {}

# Dictionary to store dependencies for each module 
# dependency_dict[fname] = dependency_set
dependeny_dict = {}

# Scan Fortran file and extract module names contained within a file
def scan_file_for_modules(file_path):
    file_path_full = os.path.join(src_list[file_path], file_path)
    with open(file_path_full, 'r') as f:
        modules = {}
        
        for line in f:
            # Look for module definitions
            module_match = module_pattern.match(line)
            if module_match:
                modules[module_match.group(1).lower()] = file_path

    return modules

# Function to scan a Fortran file and extract dependencies
def scan_fortran_file_for_deps(file_path):
    file_path_full = os.path.join(src_list[file_path], file_path)
    with open(file_path_full, 'r') as f:
        used_modules = set()
        
        for line in f:
            # Look for 'use' statements (including 'only')
            use_match = use_pattern.match(line)
            if use_match:
                used_modules.add(use_match.group(1).lower())

    return used_modules

# Walk through the current directory and its subdirectories and extract filenames and modules
for i, src_dir in enumerate(src_dirs):
    for root, _, files in os.walk(src_dir):
        for file in files:
            if file.endswith(".f90"):
                # Get the relative file path
                file_path = os.path.relpath(os.path.join(root, file), src_dir)
                file_list.append(file_path)
                obj_list[file_path] =  build_dirs[i]
                src_list[file_path] = src_dir
                module_dict.update(scan_file_for_modules(file_path))
            
# Walk through the current directory and its subdirectories and extract dependencies for each file
for file in file_list:
    dependeny_dict[file] = scan_fortran_file_for_deps(file)
            
# Write the dependencies to the .depend file
with open(DEPEND_FILE, 'w') as dep_file:
    for file in file_list:
        file_full = os.path.join(src_list[file], file)

        # Convert the source file to object file (.f90 -> .o)
        obj_file = os.path.splitext(file)[0] + ".o"
        obj_file_full = os.path.join(obj_list[file], obj_file)

        # Write the object file and its source file
        dep_file.write(f"{obj_file_full} : {file_full}")
        
        # Write the object files for the used modules if any
        if file in dependeny_dict:
            for module in dependeny_dict[file]:
                if module in module_dict.keys():
                    module_obj_file = os.path.splitext(module_dict[module])[0] + ".o"
                    module_obj_file_full = os.path.join(obj_list[module_dict[module]], module_obj_file)
                    dep_file.write(f" {module_obj_file_full}")
        dep_file.write("\n")

print(f"Dependency file {DEPEND_FILE} generated.")