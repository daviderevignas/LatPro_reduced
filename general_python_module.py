import numpy as np
import sys





















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








