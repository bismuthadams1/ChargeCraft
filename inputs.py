import pandas as pd
import numpy as np


class ReadInput:
    """
    Read and return SMILES input as list of strings
    """
    @classmethod
    def read_smiles(cls,
    smi_file: str
    ) -> list[str]:
        smiles = []
        with open(smi_file, "r") as input:
            for mols in input:
                smiles.append(mols.strip('\n'))
        
        return smiles

    