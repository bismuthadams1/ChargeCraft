
from typing import TYPE_CHECKING, ContextManager, Dict, List, Optional
from pydantic import BaseModel, Field
from contextlib import contextmanager

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
from source.storage.db import DBConformerRecordProp , DBMoleculeRecordProp
from collections import defaultdict
from sqlalchemy.orm import Session, sessionmaker


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

        mbis_dipole: Array[float] = Field(...,
        description= "mbis dipole",
        )

        mbis_quadropole: Array[float] = Field(...,
        description= "mbis quadropole",
        )

        mbis_octopole: Array[float] = Field(...,
        description= "mbis octopole",
        )

        energy: Array[float] = Field(...,
        description= "energy associated with the calculation"                             
        )

        #charge_model_charges: Dict[str, List[float]] = Field(...,
        #description= "partial charges in JSON"
        #)



        @property
        def mulliken_charges_quantity(self) -> unit.Quantity:
            return self.mulliken_charges * unit.e

        @property
        def lowdin_charges_quantity(self) -> unit.Quantity:
            return self.lowdin_charges * unit.e

        @property
        def mbis_charges_quantity(self) -> unit.Quantity:
            return self.mbis_charges * unit.e
        
        @property
        def dipole_quantity(self) -> unit.Quantity:
            return self.dipole * unit.e * unit.bohr_radius
        
        @property
        def quadropole_quantity(self) -> unit.Quantity:
            return self.quadropole * unit.e * unit.bohr_radius**2
        
        @property
        def mbis_dipole_quantity(self) -> unit.Quantity:
            return self.mbis_dipole * unit.e * unit.bohr_radius
        
        @property
        def mbis_quadropole_quantity(self) -> unit.Quantity:
            return self.mbis_quadropole * unit.e * unit.bohr_radius ** 2 

        @property
        def mbis_octopole_quantity(self) -> unit.Quantity:
            return self.mbis_quadropole * unit.e * unit.bohr_radius ** 3       
        
        @property
        def energy_quantity(self) -> unit.Quantity:
            return self.energy * unit.hartree 
        

        @classmethod
        def from_molecule(
            cls,
            molecule: "Molecule",
            conformer: unit.Quantity,
            grid_coordinates: unit.Quantity,
            esp: unit.Quantity,
            electric_field: Optional[unit.Quantity],
            esp_settings: ESPSettings,
            variables_dictionary: dict,
            energy: int
        ) -> "MoleculePropRecord":

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
            mbis_dipole = variables_dictionary["MBIS DIPOLE"]
            mbis_quadropole = variables_dictionary["MBIS QUADRUPOLE"]
            mbis_octopole= variables_dictionary["MBIS OCTOPOLE"]

            # Create the MoleculePropRecord with the additional properties
            molecule_prop_record = MoleculePropRecord(
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
                mbis_dipole= mbis_dipole,
                mbis_quadropole= mbis_quadropole,
                mbis_octopole= mbis_octopole,
                energy = energy
            )

            return molecule_prop_record

        
class MoleculePropStore(MoleculeESPStore):
    
    def __init__(self, database_path: str = "prop-store.sqlite"):
        return super().__init__(database_path = database_path)


    @classmethod
    def _db_records_to_model(cls, db_records: List[DBMoleculeRecordProp]) -> List[MoleculePropRecord]:
        """Maps a set of database records into their corresponding
        data models.

        Parameters
        ----------
        db_records
            The records to map.

        Returns
        -------
            The mapped data models.
        """
        #molecule_esp_db_records = super()._db_records_to_model(db_records: List[DBMoleculeRecord])

        # noinspection PyTypeChecker
        return [
            MoleculePropRecord(
                tagged_smiles=db_conformer.tagged_smiles,
                conformer=db_conformer.coordinates,
                grid_coordinates=db_conformer.grid,
                esp=db_conformer.esp,
                electric_field=db_conformer.field,
                esp_settings=ESPSettings(
                    basis=db_conformer.esp_settings.basis,
                    method=db_conformer.esp_settings.method,
                    grid_settings=DBGridSettings.db_to_instance(
                        db_conformer.grid_settings
                    ),
                    pcm_settings=None
                    if not db_conformer.pcm_settings
                    else DBPCMSettings.db_to_instance(db_conformer.pcm_settings)),
                mulliken_charges = db_conformer.mulliken_charges,
                lowdin_charges = db_conformer.lowdin_charges,
                mbis_charges = db_conformer.mbis_charges,
                dipole = db_conformer.dipole,
                quadropole = db_conformer.quadropole,
                mbis_dipole= db_conformer.mbis_dipole,
                mbis_quadropole= db_conformer.mbis_quadropole,
                mbis_octopole= db_conformer.mbis_octopole,
                energy = db_conformer.energy
            )
            for db_record in db_records
            for db_conformer in db_record.conformers
        ]


    @classmethod
    def _store_smiles_records(
        cls, db: Session, smiles: str, records: List[MoleculePropRecord]
    ) -> DBMoleculeRecordProp:
        """Stores a set of records which all store information for the same
        molecule.

        Parameters
        ----------
        db
            The current database session.
        smiles
            The smiles representation of the molecule.
        records
            The records to store.
        """

        existing_db_molecule = (
            db.query(DBMoleculeRecordProp).filter(DBMoleculeRecordProp.smiles == smiles).first()
        )

        if existing_db_molecule is not None:
            db_record = existing_db_molecule
        else:
            db_record = DBMoleculeRecordProp(smiles=smiles)

        # noinspection PyTypeChecker
        # noinspection PyUnresolvedReferences
        db_record.conformers.extend(
            DBConformerRecordProp(
                tagged_smiles=record.tagged_smiles,
                coordinates=record.conformer,
                grid=record.grid_coordinates,
                esp=record.esp,
                field=record.electric_field,
                grid_settings=DBGridSettings.unique(
                    db, record.esp_settings.grid_settings
                ),
                pcm_settings=None
                if not record.esp_settings.pcm_settings
                else DBPCMSettings.unique(db, record.esp_settings.pcm_settings),
                esp_settings=DBESPSettings.unique(db, record.esp_settings),
                mulliken_charges = record.mulliken_charges,
                lowdin_charges = record.lowdin_charges,
                mbis_charges = record.mbis_charges,
                dipole = record.dipole,
                quadropole = record.quadropole,
                mbis_dipole= record.mbis_dipole,
                mbis_quadropole= record.mbis_quadropole,
                mbis_octopole= record.mbis_octopole,
                energy = record.energy
            )
            for record in records
        )

        if existing_db_molecule is None:
            db.add(db_record)

        return db_record


    def store(self, *records: MoleculePropRecord):
        """Store the properties  calculated for
        a given molecule in the data store.

        Parameters
        ----------
        records
            The records to store.

        Returns
        -------
            The records as they appear in the store.
        """

        # Validate an re-partition the records by their smiles patterns.
        records_by_smiles: Dict[str, List[MoleculePropRecord]] = defaultdict(list)

        for record in records:
            record = MoleculePropRecord(**record.dict())
            smiles = self._tagged_to_canonical_smiles(record.tagged_smiles)

            records_by_smiles[smiles].append(record)

        # Store the records.
        with self._get_session() as db:
            for smiles in records_by_smiles:
                self._store_smiles_records(db, smiles, records_by_smiles[smiles])            

    def store_partial(self, tagged_smiles: str, conformer:Array, charge_model: str, charges: Array) -> None:
        ...

    def retrieve_partial(self, tagged_smiles: str, conformer:Array ) -> None:
        ...

    def retrieve(
            self,
            smiles: Optional[str] = None,
            basis: Optional[str] = None,
            method: Optional[str] = None,
            implicit_solvent: Optional[bool] = None,
        ) -> List[MoleculePropRecord]:
            """Retrieve records stored in this data store, optionally
            according to a set of filters."""

            with self._get_session() as db:
                db_records = db.query(DBMoleculeRecordProp)

                if smiles is not None:
                    smiles = self._tagged_to_canonical_smiles(smiles)
                    db_records = db_records.filter(DBMoleculeRecordProp.smiles == smiles)

                if basis is not None or method is not None or implicit_solvent is not None:
                    db_records = db_records.join(DBConformerRecordProp)

                    if basis is not None or method is not None:
                        db_records = db_records.join(
                            DBESPSettings, DBConformerRecordProp.esp_settings
                        )

                        if basis is not None:
                            db_records = db_records.filter(DBESPSettings.basis == basis)
                        if method is not None:
                            db_records = db_records.filter(DBESPSettings.method == method)

                    if implicit_solvent is not None:
                        if implicit_solvent:
                            db_records = db_records.filter(
                                DBConformerRecordProp.pcm_settings_id.isnot(None)
                            )
                        else:
                            db_records = db_records.filter(
                                DBConformerRecordProp.pcm_settings_id.is_(None)
                            )

                db_records = db_records.all()

                records = self._db_records_to_model(db_records)

                if basis:
                    records = [
                        record for record in records if record.esp_settings.basis == basis
                    ]
                if method:
                    records = [
                        record for record in records if record.esp_settings.method == method
                    ]

                return records                  

    @contextmanager
    def _get_session(self) -> ContextManager[Session]:
        session = self._session_maker()

        try:
            yield session
            session.commit()
        except BaseException as e:
            session.rollback()
            raise e
        finally:
            session.close()


    def list(self) -> List[str]:
            """Lists the molecules which exist in and may be retrieved from the
            store."""

            with self._get_session() as db:
                return [smiles for (smiles,) in db.query(DBMoleculeRecordProp.smiles).all()]