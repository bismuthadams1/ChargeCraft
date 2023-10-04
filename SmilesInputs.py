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
        """
        Reads a .smi file and returns a list of smiles strings

        Parameters
        ----------
        smi_file
            The .smi file to read.
            
        Returns
        -------
        List of smiles strings.
        """
        smiles = []
        with open(smi_file, "r") as input:
            for mols in input:
                smiles.append(mols.strip('\n'))
        
        return smiles

    