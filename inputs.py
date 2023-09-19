import pandas as pd
import numpy as np


class ReadInput:
    """
    Read and return SMILES input as list of strings
    """
    @classmethod
    def _read_smiles(cls,
    smi: str
    ) -> list[str]:
        smiles = []
        with open(smi, "r") as input:
            for mols in input:
                smiles.append(mols.strip('\n'))
        
        return smiles

    