import numpy as np
import sys
import os















##############################
#     loading functions      #
##############################


def load_hamiltonian_walks_contacts():
    tmp_hamiltonian_walks_contacts=np.loadtxt(\
        os.path.join("data/","hamiltonian_walks_contacts.csv"),delimiter=',',dtype=np.int16)
    hamiltonian_walks_contacts=[]
    for row in tmp_hamiltonian_walks_contacts:
        hamiltonian_walks_contacts.append([])
        for ii in range(len(row)//2):
            hamiltonian_walks_contacts[-1].append([row[2*ii],row[2*ii+1]])
    hamiltonian_walks_contacts=np.array(hamiltonian_walks_contacts,dtype=np.int16)
    del tmp_hamiltonian_walks_contacts
    return hamiltonian_walks_contacts

def load_hamiltonian_walks_bond_vectors():
    hamiltonian_walks_bond_vectors=np.loadtxt(\
        os.path.join("data/","hamiltonian_walks_bond_vectors.csv"),delimiter=',',dtype=np.int16)
    return hamiltonian_walks_bond_vectors

def load_MJ_tab_6():
    tmp_MJ_tab_6=np.loadtxt(os.path.join("data/","MJ_tab_6.csv"),delimiter=',')
    MJ_tab_6=np.zeros((20,20))
    for ii in range(len(tmp_MJ_tab_6)):
        MJ_tab_6[int(tmp_MJ_tab_6[ii,0])-1, int(tmp_MJ_tab_6[ii,1])-1] = tmp_MJ_tab_6[ii,2]
    del tmp_MJ_tab_6
    return MJ_tab_6



#############################
#     common variables      #
#############################

num_to_letter_map = np.array([i for i in "CMFILVWYAGTSNQDEHRKP"])
letter_to_num_map = {letter: idx for idx, letter in enumerate(num_to_letter_map)}
aa_property_map = np.array([
    1, # Cysteine (Polar Uncharged)
    0, # Methionine (Nonpolar)
    0, # Phenylalanine (Nonpolar)
    0, # Isoleucine (Nonpolar)
    0, # Leucine (Nonpolar)
    0, # Valine (Nonpolar)
    0, # Tryptophan (Nonpolar)
    1, # Tyrosine (Polar Uncharged)
    0, # Alanine (Nonpolar)
    0, # Glycine (Nonpolar) **
    1, # Threonine (Polar Uncharged)
    1, # Serine (Polar Uncharged)
    1, # Asparagine (Polar Uncharged)
    1, # Glutamine (Polar Uncharged)
    2, # Aspartic Acid (Negatively Charged)
    2, # Glutamic Acid (Negatively Charged)
    3, # Histidine (Positively Charged/Basic) **
    3, # Arginine (Positively Charged/Basic)
    3, # Lysine (Positively Charged/Basic)
    0  # Proline (Nonpolar)
])











