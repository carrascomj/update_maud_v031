"""Script to transform Maud model to the current version."""
from typing import Optional

import toml
from pydantic import BaseModel, Field

import data_model as fd


class Compartment(BaseModel):
    id: str
    name: str
    volume: float


class MetabolitePrev(BaseModel):
    id: str = Field(alias="metabolite")
    compartment: str
    balanced: bool
    name: Optional[str] = None
    inchi_key: Optional[str] = None


class ModifierPrev(BaseModel):
    modifier_type: str
    mic_id: str


class EnzymePrev(BaseModel):
    id: str
    name: str
    subunits: int = 1
    modifier: Optional[list[ModifierPrev]]


class ReactionPrev(BaseModel):
    id: str
    name: str
    reaction_mechanism: Optional[str] = "reversible_modular_rate_law"
    stoichiometry: dict[str, float]
    enzyme: list[EnzymePrev]
    water_stoichiometry: float = 0
    transported_charge: float = 0


class DrainPrev(BaseModel):
    id: str
    name: str
    stoichiometry: dict[str, float]


class ModelPrev(BaseModel):
    compartment: list[Compartment]
    reaction: list[ReactionPrev]
    metabolite: list[MetabolitePrev] = Field(alias="metabolite-in-compartment")
    drain: list[DrainPrev]

    class Config:
        allow_population_by_field_name = True


class ModelNew(BaseModel):
    compartment: list[Compartment]
    reaction: list[fd.Reaction]
    enzyme: list[fd.Enzyme]
    enzyme_reaction: list[fd.EnzymeReaction]
    metabolite: list[fd.Metabolite]
    metabolite_in_compartment: list[fd.MetaboliteInCompartment]
    allostery: list[fd.Allostery]
    competitive_inhibition: list[fd.CompetitiveInhibition]


def update_model(old_model: ModelPrev) -> ModelNew:
    """Translate a Maud old model to a new model."""
    reactions = []
    enzymes = []
    enzyme_reactions = []
    allosteries = []
    comp_inhibitions = []
    for reac in old_model.reaction:
        mechanism = (
            fd.ReactionMechanism.REVERSIBLE_MICHAELIS_MENTEN
            if reac.reaction_mechanism.startswith("reversible")
            else fd.ReactionMechanism.IRREVERSIBLE_MICHAELIS_MENTEN
        )
        reac_id = reac.id.replace("_", "")
        reactions.append(
            fd.Reaction(
                id=reac_id,
                name=reac.name,
                mechanism=mechanism,
                stoichiometry=reac.stoichiometry,
                water_stoichiometry=reac.water_stoichiometry,
                transported_charge=reac.transported_charge,
            )
        )
        for enz in reac.enzyme:
            enz_id = enz.id.replace("_", "")
            enzymes.append(fd.Enzyme(id=enz_id, name=enz.name, subunits=enz.subunits))
            enzyme_reactions.append(
                fd.EnzymeReaction(enzyme_id=enz_id, reaction_id=reac_id)
            )
            if enz.modifier is not None:
                for modifier in enz.modifier:
                    if modifier.modifier_type == "competitive_inhibitor":
                        comp_inhibitions.append(
                            fd.CompetitiveInhibition(
                                enzyme_id=enz_id,
                                reaction_id=reac_id,
                                metabolite_id=modifier.mic_id,
                                compartment_id="c",
                            )
                        )
                    else:
                        mod_type = (
                            fd.ModificationType.INHIBITION
                            if modifier.modifier_type == "allosteric_inhibitor"
                            else fd.ModificationType.ACTIVATION
                        )
                        allosteries.append(
                            fd.Allostery(
                                enzyme_id=enz_id,
                                metabolite_id=modifier.mic_id,
                                compartment="c",
                                modification_type=mod_type,
                            )
                        )
    metabolites = []
    comp_metabolites = []
    for met in old_model.metabolite:
        met_id = met.id.replace("_", "")
        metabolites.append(
            fd.Metabolite(id=met_id, name=met.name, inchi_key=met.inchi_key)
        )
        comp_metabolites.append(
            fd.MetaboliteInCompartment(
                metabolite_id=met_id,
                compartment_id=met.compartment,
                balanced=met.balanced,
            )
        )
    return ModelNew(
        compartment=old_model.compartment,
        enzyme=enzymes,
        reaction=reactions,
        enzyme_reaction=enzyme_reactions,
        metabolite=metabolites,
        metabolite_in_compartment=comp_metabolites,
        allostery=allosteries,
        competitive_inhibition=comp_inhibitions,
    )


def write_new_model(model: ModelNew, out_file: str):
    with open(out_file, "w") as f:
        toml.dump(
            model.dict(),
            f,
        )


def read_old_maud(toml_file: str):
    with open(toml_file) as f:
        data = ModelPrev.parse_obj(toml.load(f))
    return data


if __name__ == "__main__":
    data = read_old_maud("/home/georg/models/kinetic_cauto.git/trunk/data/cauto.toml")
    new_data = update_model(data)
    write_new_model(new_data, "model.toml")
