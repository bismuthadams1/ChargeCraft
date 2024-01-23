"""Add the qubekit charges and virtual sites to the database."""
import os
import sys

#sys.path.append("/Users/joshuahorton/Documents/Software/QM_ESP_Psi4")

from qubekit.molecules import Ligand, Dipole, Quadrupole
from qubekit.nonbonded.virtual_sites import VirtualSites
from chargecraft.storage.storage import MoleculePropStore
from openff.toolkit.topology import Molecule
from openff.units import unit

prop_store = MoleculePropStore(database_path="/properties_store.db")
vs = VirtualSites()

for smiles in prop_store.list():
    print(f"Running {smiles}")
    mol_props = prop_store.retrieve(smiles=smiles)
    for mol_prop in mol_props:
        off_mol = Molecule.from_mapped_smiles(mol_prop.tagged_smiles)
        off_mol.add_conformer(mol_prop.conformer_quantity)
        # convert to qubekit model
        qubekit_ligand = Ligand.from_rdkit(off_mol.to_rdkit())
        qubekit_ligand.name = smiles
        mbis_charges = mol_prop.mbis_charges
        mbis_dipoles = mol_prop.mbis_dipole
        mbis_quads = mol_prop.mbis_quadropole
        for i in range(qubekit_ligand.n_atoms):
            qube_atom = qubekit_ligand.atoms[i]
            # assign all data
            qube_atom.aim.charge = mbis_charges[i][0]
            qubekit_ligand.NonbondedForce.create_parameter(atoms=(i, ), charge=mbis_charges[i][0], epsilon=0, sigma=0)
            dipole = Dipole(x=mbis_dipoles[i][0], y=mbis_dipoles[i][1], z=mbis_dipoles[i][2])
            qube_atom.dipole = dipole
            trace = mbis_quads[i][0][0] + mbis_quads[i][1][1] + mbis_quads[i][2][2]
            trace /= 3
            # make sure we have the traceless quad tensor
            quad = Quadrupole(
                q_xx=mbis_quads[i][0][0] - trace,
                q_xy=mbis_quads[i][0][1],
                q_xz=mbis_quads[i][0][2],
                q_yy=mbis_quads[i][1][1] - trace,
                q_zz=mbis_quads[i][2][2] - trace,
                q_yz=mbis_quads[i][1][2],
            )
            qube_atom.quadrupole = quad
        qubekit_ligand = vs.run(qubekit_ligand)
        # read charges back in

        if qubekit_ligand.extra_sites.n_sites > 0:
            print(f"Molecule has sites {smiles}")
            coords = []
            charges = []
            with open(os.path.join(smiles, "xyz_with_extra_point_charges.xyz")) as xyz:
                for line in xyz.readlines()[2:]:
                    coords.extend([float(f) for f in line.split()[1:4]])
                    charges.append(float(line.split()[-1]))
            # add to db
            prop_store.store_partial(
                smiles=mol_prop.tagged_smiles,
                conformer=mol_prop.conformer,
                charges=charges ,
                charge_model="qubekit-charges"
            )
            if qubekit_ligand.extra_sites.n_sites > 0:
                # extract the vsite coords
                prop_store.store_partial(
                    smiles=mol_prop.tagged_smiles,
                    conformer=mol_prop.conformer,
                    charges=coords,
                    charge_model="qubekit-coords-Ang"
                )
        else:
            # just store back the mbis values
            prop_store.store_partial(
                smiles=mol_prop.tagged_smiles,
                conformer=mol_prop.conformer,
                charges=mbis_charges.flatten(),
                charge_model="qubekit-charges"
            )




