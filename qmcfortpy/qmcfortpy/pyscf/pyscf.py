# SPDX-FileCopyrightText: Copyright (c) 2025 Zoran Sukurma
# SPDX-License-Identifier: Apache-2.0

"""
Functions for export of the Hamiltonian from pyscf
"""

import os
import numpy as np
from scipy.linalg import sqrtm as sqrtm

from pyscf import ao2mo

def transform_h1(h1, trafo, basis, ispin1):
    """
    Transforms h1 to desired basis:
    if basis == "ao" h1 doesn't change
    """
    n = h1.shape[0]
    h1_mo = np.zeros([n,n,ispin1])
    if (basis == "ao"):
        h1_mo[:,:,0] = h1
    else:
        for spin in range(ispin1):
            h1_mo[:,:,spin] = np.einsum('pi,pq,qj->ij', trafo[:,:,spin], h1, trafo[:,:,spin])
    return h1_mo

def transform_h2(h2, trafo, basis, ispin1, ispin2):
    """
    Transforms h2 to desired basis
    if basis = "ao" h2 doesn't change
    """
    n = h2.shape[0]
    h2_mo = np.zeros([n,n,n,n,ispin2])

    if (basis == "ao"):
        h2_mo[:,:,:,:,0] = h2
    else:
        for spin in range(ispin1):
            spin1 = spin
            if (ispin1 == 1):
                spin1 = 0
            h2_mo[:,:,:,:,spin] = ao2mo.incore.full(h2, trafo[:,:,spin1])
        if (ispin2 == 3):
            h2_mo[:,:,:,:,2] = ao2mo.incore.general(h2, (trafo[:,:,0],trafo[:,:,0],trafo[:,:,1],trafo[:,:,1]))
    return h2_mo

def transform_orbitals(coeff, trafo, basis, ispin, ispin1):
    """
    Transforms orbital expansion coefficients to desired basis
    if basis = "ao" coeff doesn't change
    """
    n = coeff.shape[0]
    coeff_mo = np.zeros([n,n,ispin])
    
    if (basis == "ao"):
        coeff_mo = coeff
    else:
        for spin in range(ispin):
            spin1 = spin
            if (ispin1 == 1):
                spin1 = 0
            coeff_mo[:,:,spin] = np.einsum('pq,qr->pr', trafo[:,:,spin1], coeff[:,:,spin])
    return coeff_mo


def qmcfort_pyscf_itf(mol, method, basis="mo"):
    """
    Routine takes pyscf method object and pyscf Molecule object  to export
    Hamiltonian in given basis.
    """
    coeff = np.array(method.mo_coeff)
    n = coeff.shape[1]
    nel = mol.nelec

    if (len(coeff.shape) == 2):
        ispin = 1
        coeff = np.reshape(coeff, (n,n,ispin))
    else:
        ispin = 2
        coeff = np.moveaxis(coeff,0,-1)
     
    if (basis == "mo"):
        if (ispin == 2):
            ispin1 = 2
            ispin2 = 3
        else:
            ispin1 = 1
            ispin2 = 1
    else:
        ispin1 = 1
        ispin2 = 1
        
    h0 = [0.0]
    ovlp = mol.intor_symmetric("int1e_ovlp")
    h1 = mol.intor_symmetric("int1e_kin") + mol.intor_symmetric("int1e_nuc")
    h2 = mol.intor("int2e_sph", aosym=1)

    if (basis == "mo"):
        trafo = coeff
        trafo_orb = np.copy(trafo)
        for spin in range(ispin):
            trafo_orb[:,:,spin] = np.linalg.inv(trafo_orb[:,:,spin])
    else:
        if (basis == "ao_orth"):
            trafo = sqrtm(np.linalg.inv(ovlp))
            trafo_orb = sqrtm(ovlp)
            trafo = np.reshape(trafo,(n,n,ispin1))
            trafo_orb = np.reshape(trafo_orb,(n,n,ispin1))
        else:
            trafo = coeff
    
    h1_mo = transform_h1(h1, trafo, basis, ispin1)
    ovlp_mo = transform_h1(ovlp, trafo, basis, ispin1)
    h2_mo = transform_h2(h2, trafo, basis, ispin1, ispin2)
    coeff_mo = transform_orbitals(coeff, trafo_orb, basis, ispin, ispin1)

    try:
        F = method.get_fock()
        eigenval, co = method.eig(F, ovlp)
        eigenval = eigenval.reshape(n, ispin)
    except:
        eigenval = np.zeros([n, ispin])

    write_ham_info(n, nel, ispin, ispin1, ispin2, basis)
    np.save("ham0", np.array(h0,order="F"))
    np.save("ham1", np.array(h1_mo,order="F"))
    np.save("ham2", np.array(h2_mo,order="F"))
    np.save("overlap", np.array(ovlp_mo,order="F"))
    np.save("orbitals", np.array(coeff_mo,order="F"))
    np.save("eigenval", np.array(eigenval,order="F"))
    os.rename("ham0.npy", "ham0")
    os.rename("ham1.npy", "ham1")
    os.rename("ham2.npy", "ham2")
    os.rename("overlap.npy", "overlap")
    os.rename("orbitals.npy", "orbitals")
    os.rename("eigenval.npy", "eigenval")


def write_ham_info(n, nel, ispin, ispin1, ispin2, basis):
    """
    Writes ham_info file
    """
    print("System Info")
    print("===========")
    print("  n                     = ", str(n))
    print("  nbtot                 = ", str(n))
    print("  ispin                 = ", str(ispin))
    print("  ispin1                = ", str(ispin1))
    print("  ispin2                = ", str(ispin2))
    print("  nel                   = ", nel[0], nel[1])
    print("  basis                 = ", basis)
    print("==========================================")

    f = open("ham_info", "w")
    f.write("basis               = " + str(basis) + "\n")
    f.write("integral_mode       = " + "eri" + "\n")
    f.write("file_format         = " + "numpy" + "\n")
    f.write("n                   = " + str(n) + "\n")
    f.write("nbtot               = " + str(n) + "\n")
    f.write("ispin               = " + str(ispin) + "\n")
    f.write("ispin1              = " + str(ispin1) + "\n")
    f.write("ispin2              = " + str(ispin2) + "\n")
    #f.write("ng                  = " + str(ng) + "\n")
    f.write("nel                 = " + str(nel[0]) + "  " + str(nel[1]) + "\n")
    f.close()

    return None


def write_cas_space(casscf, tol=0.0):
    """
    Write determinants and CI coefficients of the given CAS space
    """
    n = casscf.ncas
    ne = casscf.nelecas
    cas = casscf.fcisolver.large_ci(casscf.ci, n, ne, tol, return_strs=False)
    ndets = len(cas)
    ndets_tot = np.size(casscf.ci)
    ci_arr = []
    for ci, ia, ib in cas:
        ci_arr.append(ci)
    idxs = np.argsort(np.abs(np.array(ci_arr)))
    cas_sort = []
    for idx in np.flip(idxs):
        cas_sort.append(cas[idx])
    cas = cas_sort
    print("CAS tolerance for determinants         = ", tol)
    print("CASSCF total number of determinants    = ", ndets_tot)
    print("CASSCF dets with coeff larger than tol = ", ndets)
    f = open("cas_space", "w")
    f.write("# of dets   # of active orbitals   # of active electrons"  + "\n")
    f.write("    " + str(ndets)+ "                " + str(n) + "      " + "         " + str(ne[0]) + "     " + str(ne[1]) + "\n")
    f.write("   alpha string          beta string            ci_coeff" + "\n")
    for coeff, ia, ib in cas:
        astr = "[ " + " ".join(str(i) for i in ia) + " ]"
        bstr = "[ " + " ".join(str(i) for i in ib) + " ]"
        f.write(astr + "  " + bstr + "  " + str(coeff) + "\n")
    f.close()
    return None