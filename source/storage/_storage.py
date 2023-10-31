
from typing import TYPE_CHECKING, ContextManager, Dict, List, Optional
from pydantic import BaseModel, Field
import numpy as np
from openff.recharge.esp.storage import MoleculeESPRecord, MoleculeESPStore
from openff.recharge.esp import ESPSettings
from openff.recharge.esp.storage.db import (
    DB_VERSION,
    DBBase,
    DBGridSettings,
    DBConformerRecord,
    DBESPSettings,
    DBGeneralProvenance,
    DBInformation,
    DBMoleculeRecord,
    DBPCMSettings,
    DBSoftwareProvenance,
)
from openff.recharge.esp.storage.exceptions import IncompatibleDBVersion


if TYPE_CHECKING:
    from openff.toolkit import Molecule

    Array = np.ndarray


class MoleculePropRecord(MoleculeESPRecord):
        """An extension of the MoleculeESPRecord class to store ESPs, partial charges and quadropoles. 
        """

        mulliken_charges: Array[float] = Field(...,
        description="The muliken charges associated with each atom in the conformer",
        )

        lowdin_charges: Array[float] = Field(...,
        description="The lowdin charges associated with each atom in the conformer",
        )


