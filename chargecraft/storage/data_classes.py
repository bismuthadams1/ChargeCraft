from openff.recharge.esp import ESPSettings as OldESPSettings
from pydantic import BaseModel, Field
from typing import TYPE_CHECKING, Literal, Optional, Tuple, Union
from openff.recharge.esp import DFTGridSettings
from openff.recharge.grids import GridGenerator, GridSettingsType

from pydantic import BaseModel, Field


if TYPE_CHECKING:
    from openff.toolkit import Molecule

    PositiveFloat = float
else:
    from pydantic import PositiveFloat


class DDXSettings(BaseModel):
    """A class which described the DDX model settings to include in the calculation of an ESP
    """

    solvent: Optional[str] = Field(None,
    description="The solvent to simulate. This can be 'Water' as a string",
    )

    epsilon: PositiveFloat = Field(None,
    description="the dialectric constant of the fluid")

    radii_set: Literal["uff","bondi"] = Field("uff",
    description="The type of atomic radii to use when computing the molecular "
        "cavity.",
        )

    ddx_model: Literal["PCM","COSMO","LPB"] = Field("PCM",
    description="The available solvation models",
    )
    
class PCMSettings(BaseModel):
    """A class which describes the polarizable continuum model (PCM)
    to include in the calculation of an ESP.
    """

    solver: Literal["CPCM", "IEFPCM"] = Field("CPCM", description="The solver to use.")

    solvent: str = Field("Water",
     description="The solvent to simulate. This controls the dielectric constant "
        "of the model.",
    )

    radii_model: Literal["Bondi", "UFF", "Allinger"] = Field(
        "Bondi",
        description="The type of atomic radii to use when computing the molecular "
        "cavity.",
    )
    radii_scaling: bool = Field(
        True, description="Whether to scale the atomic radii by a factor of 1.2."
    )

    cavity_area: PositiveFloat = Field(
        0.3, description="The average area of the surface partition for the cavity."
    )

class ESPSettings(BaseModel):
    """A class which contains the settings to use in an ESP calculation."""

    basis: str = Field(
        "6-31g*", description="The basis set to use in the ESP calculation."
    )
    method: str = Field("hf", description="The method to use in the ESP calculation.")

    grid_settings: GridSettingsType = Field(
        ...,
        description="The settings to use when generating the grid to generate the "
        "electrostatic potential on.",
    )

    pcm_settings: Optional[PCMSettings] = Field(
        None,
        description="The settings to use if including a polarizable continuum "
        "model in the ESP calculation.",
    )

    psi4_dft_grid_settings: DFTGridSettings = Field(
        DFTGridSettings.Default,
        description="The DFT grid settings to use when performing computations with "
        "Psi4.",
    )

    ddx_settings: Optional[DDXSettings] = Field(None,
    description="The settings to use if including the DDX model"
        "model in the ESP calculation.",
    )


