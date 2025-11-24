import os
import numpy as np

def clean_files():
    if os.path.exists("ham0"):
        os.remove("ham0")
    if os.path.exists("ham1"):
        os.remove("ham1")
    if os.path.exists("ham2"):
        os.remove("ham2")
    if os.path.exists("orbitals"):
        os.remove("orbitals")
    if os.path.exists("overlap"):
        os.remove("overlap")
    if os.path.exists("eigenval"):
        os.remove("eigenval")
    if os.path.exists("ham_info"):
        os.remove("ham_info")
    if os.path.exists("cas_space"):
        os.remove("cas_space")
    if os.path.exists("qmcfort_pos"):
        os.remove("qmcfort_pos")
    return None