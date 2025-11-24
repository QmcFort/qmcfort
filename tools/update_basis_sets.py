import basis_set_exchange as bse
import sys
import os
import shutil

if (len(sys.argv) > 1):
    dirname = sys.argv[1]
else:
    dirname = "data/basis_sets"

basis_sets = ["STO-3G", "STO-6G", "3-21G", "6-31G", "6-311G", "6-31G*", "6-311G*", "6-31G**", "6-311G**",
              "cc-pVDZ", "cc-pVTZ", "cc-pVQZ", "cc-pV5Z",
              "aug-cc-pVDZ", "aug-cc-pVTZ", "aug-cc-pVQZ", "aug-cc-pV5Z",
              "ano-pVDZ", "ano-pVTZ", "ano-pVQZ", "ano-pV5Z",
              "aug-ano-pVDZ", "aug-ano-pVTZ", "aug-ano-pVQZ", "aug-ano-pV5Z",
              "def2-SVP", "def2-TZVP", "def2-QZVP"]

atoms = ["H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl", "Ar"]

if (os.path.exists(dirname) and os.path.isdir(dirname)):
    shutil.rmtree(dirname)
os.makedirs(dirname)

for basis_set in basis_sets:
    os.mkdir(dirname + "/" + basis_set)
    for atom in atoms:
        fname = dirname + "/" + basis_set + "/" + atom
        file = open(fname, "w")
        file.write(bse.get_basis(basis_set, elements=atom, fmt="Gaussian94"))
        file.close()
