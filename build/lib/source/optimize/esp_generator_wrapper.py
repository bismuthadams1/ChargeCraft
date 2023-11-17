from source.optimize.openff_psi4_gen import Psi4Generate
from source.storage.storage import MoleculePropRecord, MoleculePropStore
from openff.toolkit import Molecule
from openff.recharge.esp import ESPSettings
from openff.recharge.esp.psi4 import Psi4ESPGenerator
from openff.recharge.esp.exceptions import Psi4Error
from tqdm import tqdm
from openff.recharge.esp.storage import MoleculeESPRecord, MoleculeESPStore
from rdkit.Chem import rdmolfiles
import qcengine
from openff.recharge.grids import GridSettingsType, GridGenerator
from openff.units import unit
from qcelemental.models import Molecule as QCMolecule
from qcelemental.models.common_models import Model
from qcelemental.models.procedures import OptimizationInput, QCInputSpecification
import copy
import numpy as np


class ESPGenerator:
    """"
    Class for Generating ESPs which wraps around the Psi4ESPGenerator class and handles errors.   
    """

    def __init__(
              self, 
              molecule: "Molecule",
              conformers: list[unit.Quantity],
              esp_settings: "ESPSettings",
              grid_settings: "GridSettingsType",
              ncores: int | None = None,
              memory: int | None = None
              ) -> None:
        self.conformers = conformers
        self.molecule = molecule
        #self.rd_molecule = self.molecule.to_rdkit()
        self.esp_settings = esp_settings
        self.grid_settings = grid_settings
        self.qc_data_store = MoleculeESPStore()
        self.records = []
        self.ncores = ncores
        self.memory = memory

    @property
    def ncores(self):
        return self._ncores
    
    @ncores.setter
    def ncores(self, value):
        if value is None:
            self._ncores = qcengine.get_config().ncores
        else:
            self._ncores = value
    
    @property
    def memory(self):
        return self._memory
    
    @memory.setter
    def memory(self, value):
        if value is None:
            self._memory =  qcengine.get_config().memory# * 0.9
        else:
            self._memory = value


    def run_esps(self) -> None:
        """
        Run psi4 to generate the ESPs, the function loops through the conformers and handles errors.
        Appends the outputs to the sqlfile.

        Parameters
        ----------

        Returns
        -------
            The contents of the input file.
        """

        for conf_no in range(self.molecule.n_conformers):
            print(f'conformer {conf_no} for {self.molecule.to_smiles()}')
            qc_mol = self.molecule.to_qcschema(conformer=conf_no)
            #run a ff optimize for each conformer to make sure the starting structure is sensible
            xtb_opt_mol = self._xtb_ff_opt(qc_mol)
            # run the psi4 HF opt starting from the final conformer
            hf_opt_mol = self._psi4_opt(qc_mol=xtb_opt_mol)
            opt_molecule = Molecule.from_qcschema(hf_opt_mol)
        
            conformer, grid, esp, electric_field, E = Psi4ESPGenerator.generate(
               molecule=opt_molecule, conformer=opt_molecule.conformers[0], settings=self.esp_settings, minimize=False
            )
           
            record = MoleculeESPRecord.from_molecule(
                    self.molecule, conformer, grid, esp, electric_field, self.esp_settings, E
                )
            # push records to the db
            self.records.append(record)
        self.qc_data_store.store(*self.records)

    def fetch_data(self) -> list[MoleculeESPRecord]:
        """
        Fetch the data from the qc_data_store. 

        Parameters
        ----------

        Returns
        -------
            The ESP contents 
        """
        self.qc_data_store.store(*self.records)

        # Retrieve the stored properties.
        return self.qc_data_store.retrieve()
    
    def _xtb_ff_opt(self, qc_mol: QCMolecule) -> QCMolecule:
        """
        Runs an xtb ff optimisation on the conformer using the qc_engine wrapper. 

        Parameters
        ----------
        conformer_no
            The conformer number in the molecule to run the optimisation on. 

        Returns
        -------
            The ff optimised conformer. 
        """
        xtb_model = Model(method="gfn2-xtb", basis=None)
        keywords = {"verbosity": "muted"}
        return self._qcengine_opt(
            qc_mol=qc_mol, model=xtb_model, program="xtb", spec_keywords=keywords
        )

    def _qcengine_opt(self, qc_mol: QCMolecule, model: Model, program: str, spec_keywords: dict[str, str]) -> QCMolecule:
        """
        A general function to run an optimisation via qcengine.
        """
        spec = QCInputSpecification(model=model, keywords=spec_keywords, driver="gradient")
        opt_spec = OptimizationInput(
            initial_molecule=qc_mol,
            input_specification=spec,
            keywords={"coordsys": "dlc", 
                      "program": program
                      }
                    
        )
        opt = qcengine.compute_procedure(opt_spec, "geometric", local_options = { "memory": self.memory, "ncores": self.ncores})
        print(opt)
        return opt.final_molecule

    def _psi4_opt(self, qc_mol: QCMolecule) -> QCMolecule:
        """
        Run an geometry optimivsation using PSI4 and geometric.

        Parameters
        ----------
        qc_mol:
            The qcelemental version of the molecule to be optimised.

        Returns
        -------
            The optimised qcelemental molecule
        """
        hf_model = Model(method="hf", basis="6-31G*")
        # do want anything here like density fitting?
        keywords = {}
        return self._qcengine_opt(
            qc_mol=qc_mol, model=hf_model, program="psi4", spec_keywords=keywords
        )


    def _generate_grid(self, 
                       conformer: unit.Quantity) -> unit.Quantity:
        """
        Generates the grid for the ESP. 

        Parameters
        ----------
        conformer
            The conformer which the grid needs to be generated on.

        Returns
        -------
            The grid. 
        """
        grid = GridGenerator.generate(self.molecule, conformer, self.grid_settings)
        #returns grid in Angstrom
        return grid


class PropGenerator(ESPGenerator):
    """"
    Class for Generating ESPs and properties which wraps around the Psi4Generate class and handles errors.   
    """
    def __init__(self,
                 molecule: "Molecule",
                 conformers: list[unit.Quantity],
                 esp_settings: "ESPSettings",
                 grid_settings: "GridSettingsType",
                 ncores: int | None = None,
                 memory: int | None = None,
                 prop_data_store = MoleculePropStore()
                 ) -> None:
        super().__init__(
                   molecule,
                   conformers,
                   esp_settings,
                   grid_settings,
                   ncores,
                   memory)
        self.prop_data_store = prop_data_store

        
    
    def run_props(self) -> None:
        """
        Run psi4 to generate the ESPs and the properties, the function loops through the conformers and handles errors.
        Appends the outputs to the sqlfile.

        Parameters
        ----------

        Returns
        -------
            The contents of the input file.
        """
        for conf_no in range(self.molecule.n_conformers):
            print(f'conformer {conf_no} for {self.molecule.to_smiles()}')
            # #The default dynamic level is 1, we've made it higher to
            dynamic_level = 5
            qc_mol = self.molecule.to_qcschema(conformer=conf_no)
            #run a ff optimize for each conformer to make sure the starting structure is sensible
            xtb_opt_mol = self._xtb_ff_opt(qc_mol)
            #hf QM optimisation           
            hf_opt_mol = self._psi4_opt(qc_mol=xtb_opt_mol)
            #ensure molecule is not orientated
            qc_mol_opt = QCMolecule(**hf_opt_mol.dict(exclude={"fix_com", "fix_orientation"}), fix_com=True, fix_orientation=True)
            #QC geometries a given as bohr in the qcelement .geometry attribute, conversion to angstrom is handled in the GridGenerator.generate() method. 
            grid = self._generate_grid((qc_mol_opt.geometry * unit.bohr))
            try:
                 conformer, grid, esp, electric_field, variables_dictionary, E  = self._prop_generator_wrapper(conformer = qc_mol_opt, dynamic_level = dynamic_level, grid = grid)
            except Psi4Error:
                 #if this conformer after a few attempts (contained in _esp_generator_wrapper function) the move to the next conformer.
                 continue
            
            record = MoleculePropRecord.from_molecule(
                    self.molecule, conformer, grid, esp, electric_field, self.esp_settings, variables_dictionary, E
            )
            print(*record)
            self.records.append(record)
            print(variables_dictionary)
        self.prop_data_store.store(*self.records)

    def _prop_generator_wrapper(self, 
                                conformer: "QCMolecule", 
                                dynamic_level: int, 
                                grid: unit.Quantity,
                                error_level: int = 0) -> tuple[unit.Quantity, unit.Quantity, unit.Quantity, unit.Quantity] | Psi4Error:
        """
        Wrapper around the Psi4ESP generator which slowly increases the dynamic level if the calculation has problems

        Parameters
        ----------
        conformer
            The conformer of the molecule to generate the ESP for.
        dynamic_level
            This corresponds to the DYNAMIC_LEVEL keyword in the psi4 calculation. 
        grid
            ESP gridpoints.

        Returns
        -------
            Conformer, grid, esp, and electric field OR a Psi4Error
            
        """
        # run through different error options, slowly escalate.
        try:
            xyz, grid, esp, electric_field, variables_dictionary, E = Psi4Generate.get_properties(
                            molecule = self.molecule,
                            conformer = conformer,
                            grid = grid,
                            settings = self.esp_settings,
                            dynamic_level = dynamic_level
                    )
            return xyz, grid, esp, electric_field, variables_dictionary, E
        #Error handling, this can probably be developed. There shouldn't be any issues since the geometry will have already be optmized with geometric. This can be kept for future error handling design
        except Psi4Error:
            if error_level == 0:
                error_level += 1
                dynamic_level += 1
                grid, esp, electric_field, variables_dictionary, E = self._prop_generator_wrapper(conformer, dynamic_level, error_level)
                return xyz, grid, esp, electric_field, variables_dictionary, E
            elif error_level == 1:
                error_level += 1
                dynamic_level += 1
                grid, esp, electric_field, variables_dictionary, E = self._prop_generator_wrapper(conformer, dynamic_level, error_level)
                return xyz, grid, esp, electric_field, variables_dictionary, E
            else:
                raise Psi4Error 
