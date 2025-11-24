# pyscf_itf.py

from pyscf import gto, scf
from qmcfortpy.pyscf import qmcfort_pyscf_itf

mol = gto.M(atom = "qmcfort.xyz",
            basis = "cc-pvdz",
            spin = 0)

mf = scf.RHF(mol)
mf.kernel()

qmcfort_pyscf_itf(mol, mf, basis="ao_orth")

