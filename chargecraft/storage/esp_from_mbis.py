import numpy as np
from openff.units import unit

AU_ESP = unit.atomic_unit_of_energy / unit.elementary_charge

class ESPCalculator:
    def __init__(self):
        self.ke = 1 / (4 * np.pi * unit.epsilon_0)  # 1/vacuum_permittivity, 1/(e**2 * a0 *Eh)

    def assign_esp(self, monopoles: np.ndarray, 
                   dipoles: np.ndarray, 
                   quadropules: np.ndarray, 
                   grid: unit.Quantity, 
                   coordinates: unit.Quantity) -> list:
        """Assign charges according to charge model selected

        Parameters
        ----------
        monopoles : np.ndarray
            Array of monopoles.
        dipoles : np.ndarray
            Array of dipoles.
        quadropules : np.ndarray
            Array of quadrupoles.
        grid : unit.Quantity
            Grid coordinates.
        coordinates : unit.Quantity
            Atom coordinates.

        Returns
        -------
        list of partial charges 
        """
        monopoles_quantity = monopoles.to(unit.e)
        dipoles_quantity = dipoles.to(unit.e*unit.angstrom)
        quadropoles_quantity = quadropules.to(unit.e*unit.angstrom*unit.angstrom)
        coordinates_ang = coordinates.to(unit.angstrom)

        monopole_esp = self._calculate_esp_monopole_au(grid_coordinates=grid, atom_coordinates=coordinates_ang, charges=monopoles_quantity)
        dipole_esp = self._calculate_esp_dipole_au(grid_coordinates=grid, atom_coordinates=coordinates_ang, dipoles=dipoles_quantity)
        quadrupole_esp = self._calculate_esp_quadropole_au(grid_coordinates=grid, atom_coordinates=coordinates_ang, quadrupoles=quadropoles_quantity)

        return (monopole_esp + dipole_esp + quadrupole_esp).m.flatten().tolist(), grid.m.tolist()

    def _calculate_esp_monopole_au(self, 
                                  grid_coordinates: unit.Quantity, 
                                  atom_coordinates: unit.Quantity, 
                                  charges: unit.Quantity) -> unit.Quantity:
        """Generate the esp from the on atom monopole

        Parameters
        ----------
        grid_coordinates : unit.Quantity
            Grid on which to build the esp on.
        atom_coordinates : unit.Quantity
            Coordinates of atoms to build the esp.
        charges : unit.Quantity
            Monopole or charges.

        Returns
        -------
        monopole_esp : unit.Quantity
            Monopole esp.
        """
        charges = charges.flatten()
        grid_coordinates = grid_coordinates.reshape((-1, 3)).to(unit.bohr)
        atom_coordinates = atom_coordinates.reshape((-1, 3)).to(unit.bohr)

        displacement = grid_coordinates[:, None, :] - atom_coordinates[None, :, :]
        distance = np.linalg.norm(displacement.m, axis=-1) * unit.bohr
        inv_distance = 1 / distance

        esp = self.ke * np.sum(inv_distance * charges[None, :], axis=1)
        return esp.to(AU_ESP)

    def _calculate_esp_dipole_au(self, 
                                grid_coordinates: unit.Quantity, 
                                atom_coordinates: unit.Quantity, 
                                dipoles: unit.Quantity) -> unit.Quantity:
        """Generate the esp from the on atom dipoles

        Parameters
        ----------
        grid_coordinates : unit.Quantity
            Grid on which to build the esp on.
        atom_coordinates : unit.Quantity
            Coordinates of atoms to build the esp.
        dipoles : unit.Quantity
            Dipoles or charges.

        Returns
        -------
        dipoles_esp : unit.Quantity
            Monopole esp.
        """
        dipoles = dipoles.to(unit.e * unit.bohr)
        if dipoles.ndim == 1:
            dipoles = dipoles.reshape((-1, 3))

        grid_coordinates = grid_coordinates.reshape((-1, 3)).to(unit.bohr)
        atom_coordinates = atom_coordinates.reshape((-1, 3)).to(unit.bohr)

        displacement = grid_coordinates[:, None, :] - atom_coordinates[None, :, :]
        distance = np.linalg.norm(displacement.m, axis=-1) * unit.bohr
        inv_distance_cubed = 1 / distance**3
        dipole_dot = np.sum(displacement * dipoles[None, :, :], axis=-1)

        esp = self.ke * np.sum(inv_distance_cubed * dipole_dot, axis=1)
        return esp.to(AU_ESP)

    def _calculate_esp_quadropole_au(self, 
                                    grid_coordinates: unit.Quantity, 
                                    atom_coordinates: unit.Quantity, 
                                    quadrupoles: unit.Quantity) -> unit.Quantity:
        """Generate the esp from the on atom quadrupoles

        Parameters
        ----------
        grid_coordinates : unit.Quantity
            Grid on which to build the esp on.
        atom_coordinates : unit.Quantity
            Coordinates of atoms to build the esp.
        quadrupoles : unit.Quantity
            Quadrupoles or charges.

        Returns
        -------
        quadrupoles_esp : unit.Quantity
            Monopole esp.
        """
        if quadrupoles.ndim == 1:
            quadrupoles = quadrupoles.reshape((-1, 3, 3))
        quadrupoles = self._detrace(quadrupoles.to(unit.e * unit.bohr * unit.bohr))
        grid_coordinates = grid_coordinates.reshape((-1, 3)).to(unit.bohr)
        atom_coordinates = atom_coordinates.reshape((-1, 3)).to(unit.bohr)

        displacement = grid_coordinates[:, None, :] - atom_coordinates[None, :, :]
        distance = np.linalg.norm(displacement.m, axis=-1) * unit.bohr
        inv_distance = 1 / distance

        quadrupole_dot_1 = np.sum(quadrupoles[None, :, :] * displacement[:, :, None], axis=-1)
        quadrupole_dot_2 = np.sum(quadrupole_dot_1 * displacement, axis=-1)
        esp = self.ke * np.sum((3 * quadrupole_dot_2 * (1 / 2 * inv_distance**5)), axis=-1)

        return esp.to(AU_ESP)

    def _detrace(self, quadrupoles: unit.Quantity) -> unit.Quantity:
        """Make sure we have the traceless quadrupole tensor.

        Parameters
        ----------
        quadrupoles : unit.Quantity
            Quadrupoles.

        Returns
        -------
        unit.Quantity
            Detraced quadrupoles.
        """
        quadrupoles = quadrupoles.m
        for i in range(quadrupoles.shape[0]):
            trace = np.trace(quadrupoles[i])
            trace /= 3
            quadrupoles[i][0][0] -= trace
            quadrupoles[i][1][1] -= trace
            quadrupoles[i][2][2] -= trace

        return quadrupoles * unit.e * unit.bohr * unit.bohr
    
#TODO: Finalize Electric Field method
    
        # def esp_to_electric_field(self, esp: np.ndarray, grid: unit.Quantity) -> np.ndarray:
        # """
        # Convert ESP to electric field by computing the gradient.

        # Parameters
        # ----------
        # esp : np.ndarray
        #     Electrostatic potential values.
        # grid : unit.Quantity
        #     Grid coordinates.

        # Returns
        # -------
        # electric_field : np.ndarray
        #     Electric field values.
        # """
        # grid = grid.m  # Convert to base units without units for gradient calculation
        # grid_spacing = np.array([np.diff(grid[:, i]).mean() for i in range(3)])  # Assuming uniform grid spacing
        # grad_x, grad_y, grad_z = np.gradient(esp.reshape(grid.shape[:-1]), *grid_spacing)
        # electric_field = -np.array([grad_x, grad_y, grad_z])
        # return electric_field * (unit.hartree / unit.e / unit.bohr)  # Convert to appropriate units
