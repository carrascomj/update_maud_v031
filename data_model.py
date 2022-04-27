from enum import Enum
from typing import Dict, List, Optional, Union

from pydantic import BaseModel
from pydantic import Field as field
from pydantic import validator

ID_SEPARATOR = "_"


class ReactionMechanism(str, Enum):
    """Possible reaction mechanisms."""

    REVERSIBLE_MICHAELIS_MENTEN = 1
    IRREVERSIBLE_MICHAELIS_MENTEN = 2
    DRAIN = 3


class ModificationType(str, Enum):
    """Possible modification types."""

    ACTIVATION = 1
    INHIBITION = 2


class KMConfig:
    """Config allowing the KineticModel class to contain pandas objects."""

    arbitrary_types_allowed = True


class Metabolite(BaseModel):
    """Maud representation of a metabolite."""

    id: str
    name: Optional[str]
    inchi_key: Optional[str]

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v


class Enzyme(BaseModel):
    """Maud representation of an enzyme."""

    id: str
    name: Optional[str]
    subunits: int

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v

    @validator("subunits")
    def subunits_must_be_positive(cls, v):
        """Check that the subunits attribute is biologically possible."""
        assert v > 0
        return v


class Compartment(BaseModel):
    """Maud representation of an intra-cellular compartment.

    For example, cytosol or mitochondria.

    """

    id: str
    name: Optional[str]
    volume: float

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v


class Reaction(BaseModel):
    """Maud representation of a chemical reaction."""

    id: str
    name: Optional[str]
    mechanism: ReactionMechanism
    stoichiometry: Dict[str, float]
    water_stoichiometry: float
    transported_charge: float

    @validator("id")
    def id_must_not_contain_seps(cls, v):
        """Check that the id doesn't contain ID_SEPARATOR."""
        assert ID_SEPARATOR not in v
        return v

    @validator("stoichiometry")
    def stoichiometry_must_be_non_zero(cls, v):
        """Check that the stoichiometry is not zero."""
        assert v != 0
        return v


class MetaboliteInCompartment(BaseModel):
    """Maud representation of a metabolite/compartment pair.

    This is needed because metabolites often exist in multiple compartments, and
    the concentration in each one is important.

    A metabolite may also be "balanced" (i.e. not be consumed or produced at
    steady state) in one compartment but not in another.

    """

    metabolite_id: str
    compartment_id: str
    balanced: bool


class EnzymeReaction(BaseModel):
    """Maud representation of an enzyme/reaction pair.

    This is needed because some enzymes catalyse multiple reactions.

    """

    enzyme_id: str
    reaction_id: str


class Allostery(BaseModel):
    """Maud representation of an allosteric modification."""

    enzyme_id: str
    metabolite_id: str
    compartment_id: str
    modification_type: ModificationType


class CompetitiveInhibition(BaseModel):
    """Maud representation of a competitive inhibition."""

    enzyme_id: str
    reaction_id: str
    metabolite_id: str
    compartment_id: str


class Phosphorylation(BaseModel):
    """Maud representation of a phosphorylation modification."""

    enzyme_id: str
    modification_type: ModificationType
