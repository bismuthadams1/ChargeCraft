
from typing import TYPE_CHECKING, ContextManager, Dict, List, Optional
from pydantic import BaseModel, Field

import numpy as np
from openff.units import unit
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
else:
    from openff.recharge.utilities.pydantic import Array, wrapped_array_validator


class MoleculePropRecord(MoleculeESPRecord):
        """An extension of the MoleculeESPRecord class to store ESPs, partial charges and quadropoles. 
        """

        #create a series of statric variables
        mulliken_charges: Array[float] = Field(...,
        description="The muliken charges associated with each atom in the conformer",
        )

        lowdin_charges: Array[float] = Field(...,
        description="The lowdin charges associated with each atom in the conformer",
        )

        mbis_charges: Array[float] = Field(...,
        description="The lowdin charges associated with each atom in the conformer",
        )

        dipole: Array[float] = Field(...,
        description="The molecular dipole",
        )

        quadropole: Array[float] = Field(...,
        description= "molecular quadropole",
        )

        @property
        def mulliken_charges(self) -> unit.Quantity:
            return self.mulliken_charges * unit.e

        @property
        def lowdin_charges(self) -> unit.Quantity:
            return self.lowdin_charges * unit.e

        @property
        def mbis_charges(self) -> unit.Quantity:
            return self.mbis_charges * unit.e
        
        @property
        def dipole(self) -> unit.Quantity:
            return self.dipole * unit.e * unit.bohr_radius
        
        @property
        def quadropole(self) -> unit.Quantity:
            return self.quadropole * unit.e * unit.bohr_radius**2
        
        @classmethod
        def from_molecule(
            cls,
            molecule: "Molecule",
            conformer: unit.Quantity,
            grid_coordinates: unit.Quantity,
            esp: unit.Quantity,
            electric_field: Optional[unit.Quantity],
            esp_settings: ESPSettings,
            variables_dictionary: dict
        ) -> "MoleculeESPRecord":

            """Creates a new ``MoleculeESPRecord`` from an existing molecule
            object, taking care of creating the InChI and SMARTS representations.

            Parameters
            ----------
            molecule
                The molecule to store in the record.
            conformer
                The coordinates [Angstrom] of this conformer with shape=(n_atoms, 3).
            grid_coordinates
                The grid coordinates [Angstrom] which the ESP was calculated on
                with shape=(n_grid_points, 3).
            esp
                The value of the ESP [Hartree / e] at each of the grid coordinates
                with shape=(n_grid_points, 1).
            electric_field
                The value of the electric field [Hartree / (e . a0)] at each of
                the grid coordinates with shape=(n_grid_points, 3).
            esp_settings
                The settings used to generate the ESP stored in this record.
            variables_dictionary
                
            Returns
            -------
                The created record.
            """

            # Call the parent class's from_molecule method using super()
            molecule_esp_record = super().from_molecule(
                molecule, conformer, grid_coordinates, esp, electric_field, esp_settings
            )

            # Unpack the variables_dictionary and add them to the molecule prop record
            mulliken_charges = variables_dictionary["MULLIKEN_CHARGES"]
            lowdin_charges = variables_dictionary["LOWDIN_CHARGES"]
            mbis_charges = variables_dictionary["MBIS CHARGES"]
            dipole = variables_dictionary["HF DIPOLE"]
            quadropole = variables_dictionary["HF QUADRUPOLE"]

            # Create the MoleculePropRecord with the additional properties
            molecule_prop_record = cls(
                tagged_smiles=molecule_esp_record.tagged_smiles,
                conformer=molecule_esp_record.conformer,
                grid_coordinates=molecule_esp_record.grid_coordinates,
                esp=molecule_esp_record.esp,
                electric_field=molecule_esp_record.electric_field,
                esp_settings=molecule_esp_record.esp_settings,
                mulliken_charges=mulliken_charges,
                lowdin_charges=lowdin_charges,
                mbis_charges=mbis_charges,
                dipole=dipole,
                quadropole=quadropole,
            )

            return molecule_prop_record






        
class MoleculePropertiesStore(MoleculeESPStore):
    
    @classmethod
    def _db_records_to_model(cls, db_records: List[DBMoleculeRecord]) -> List[MoleculeESPRecord]:
         return super()._db_records_to_model(db_records)