from openff.recharge.esp import ESPSettings as OldESPSettings
from pydantic import BaseModel, Field
from typing import TYPE_CHECKING, Literal, Optional, Tuple, Union
from openff.recharge.esp import PCMSettings as OldPCMSettings

from pydantic import BaseModel, Field


if TYPE_CHECKING:
    from openff.toolkit import Molecule

    PositiveFloat = float
else:
    from pydantic import PositiveFloat


class DDXSettings(BaseModel):
    """A class which described the DDX model settings to include in the calculation of an ESP
    """

    solvent: Union[str,PositiveFloat] = Field("Water",
    description="The solvent to simulate. This can be 'Water' as a string or a positive float "
                "representing a custom dielectric constant.",
    )

    radii_set: Literal["uff","bondi"] = Field("uff",
    description="The type of atomic radii to use when computing the molecular "
        "cavity.",
        )

    ddx_model: Literal["PCM","COSMO","LPB"] = Field("PCM",
    description="The available solvation models",
    )

class ESPSettings(OldESPSettings):

    ddx_settings: Optional[DDXSettings] = Field(None,
    description="The settings to use if including the DDX model"
        "model in the ESP calculation.",
    )

class PCMSettings(OldPCMSettings):

    solvent: str = Field("Water",
     description="The solvent to simulate. This controls the dielectric constant "
        "of the model.",
    )

