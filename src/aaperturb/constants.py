"""constants.py
Contains constants used thoughout the modules.
"""

#This not used anywhere yet
BOND_TYPES_RDKIT = {"UNSPECIFIED": 0,
                    "SINGLE": 1,
                    "DOUBLE": 2,
                    "TRIPLE": 3,
                    "QUADRUPLE": 4,
                    "QUINTUPLE": 5,
                    "HEXTUPLE": 6,
                    "ONEANDAHALF": 7,
                    "TWOANDAHALF": 8,
                    "THREEANDAHALF": 9,
                    "FOURANDAHALF": 10,
                    "FIVEANDAHALF": 11,
                    "AROMATIC": 12,
                    "IONIC": 13,
                    "HYDROGEN": 14,
                    "THREECENTER": 15,
                    "DATIVEONE": 16,
                    "DATIVE": 17,
                    "DATIVEL": 18,
                    "DATIVER": 19,
                    "OTHER": 20,
                    "ZERO": 21}

#amino acid 1 letter codes to 3
restype_1to3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "B": "ASX",    # asp/asn
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    #"U": "SEC",    # Selenocysteine, no need, uncomment if arises
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
    "Z": "GLX"       #glu/gln
}


#ATOMS OF AMINO ACIDS
ALA = ["N", "CA",  "CB"]
ARG = ["N", "CA",  "CB", "CG", "CD", "NE", "CZ", "NH1", "NH2"]
ASN = ["N", "CA",  "CB", "CG", "OD1", "ND2"]
ASP = ["N", "CA",  "CB", "CG", "OD1", "OD2"]
CYS = ["N", "CA",  "CB", "SG"]
GLU = ["N", "CA",  "CB", "CG", "CD", "OE1", "OE2"]
GLN = ["N", "CA",  "CB", "CG", "CD", "OE1", "NE2"]
GLY = [""]
HIS = ["N", "CA",  "CB", "CG", "ND1", "CE1", "NE2", "CD2"]
ILE = ["N", "CA",  "CB", "CG1", "CG2", "CD1"]
LEU = ["N", "CA",  "CB", "CG", "CD1", "CD2"]
LYS = ["N", "CA",  "CB", "CG", "CD", "CE", "NZ"]
MET = ["N", "CA",  "CB", "CG", "SD", "CE"]
PHE = ["N", "CA",  "CB", "CG", "CD1", "CE1", "CZ", "CE2", "CD2"]
PRO = ["N", "CA",  "CB", "CG", "CD"]
SER = ["N", "CA",  "CB", "OG"]
THR = ["N", "CA",  "CB", "OG1", "CG2"]
TRP = ["N", "CA",  "CB", "CG", "CD1", "NE1", "CE2", "CD2", "CE3", "CZ3", "CH2", "CZ2"]
TYR = ["N", "CA",  "CB", "CD1", "CE1", "CZ", "CE2", "CD2", "OH"]
VAL = ["N", "CA",  "CB", "CG1", "CG2"]

#Convinient
AMINO_LST = {
                "ALA": ALA, "ARG": ARG, "ASN": ASN, "ASP": ASP, "CYS": CYS,
                "GLU": GLU, "GLN": GLN, "GLY": GLY, "HIS": HIS, "ILE": ILE,
                "LEU": LEU, "LYS": LYS, "MET": MET, "PHE": PHE, "PRO": PRO,
                "SER": SER, "THR": THR, "TRP": TRP, "TYR": TYR, "VAL": VAL
}

#for each amino acid sidechain there are some predefined axes for torsions.
ROT_AXES = {
                "ALA": [("CA", "CB")],
                "ARG": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "NE"), ("NE", "CZ")],
                "ASN": [("CA", "CB"), ("CB", "CG")],
                "ASP": [("CA", "CB"), ("CB", "CG")],
                "CYS": [("CA", "CB")],
                "GLU": [("CA", "CB"), ("CB", "CG"), ("CG", "CD")],
                "GLN": [("CA", "CB")],
                #"GLY : [("CA", "CB")],
                "HIS": [("CA", "CB"), ("CB", "CG")],
                "ILE": [("CA", "CB"), ("CB", "CG1")],
                "LEU": [("CA", "CB"), ("CB", "CG")],
                "LYS": [("CA", "CB"), ("CB", "CG"), ("CG", "CD"), ("CD", "CE")],
                "MET":[("CA", "CB"), ("CB", "CG"), ("CG", "CD")], 
                "PHE": [("CA", "CB"), ("CB", "CG")],
                #"PRO" : [("CA", "CB")],
                "SER": [("CA", "CB")],
                "THR": [("CA", "CB")],
                "TRP":[("CA", "CB"), ("CB", "CG")], 
                "TYR": [("CA", "CB"), ("CB", "CG")],
                "VAL":[("CA", "CB")]}