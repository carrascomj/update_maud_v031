"""Script to transform Maud model to the current version."""
import json
from typing import Optional

import click
import maud.data_model.kinetic_model as md
import toml
from pydantic import BaseModel, Field


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
    reaction: list[md.Reaction]
    enzyme: list[md.Enzyme]
    enzyme_reaction: list[md.EnzymeReaction]
    metabolite: list[md.Metabolite]
    metabolite_in_compartment: list[md.MetaboliteInCompartment]
    allostery: list[md.Allostery]
    competitive_inhibition: list[md.CompetitiveInhibition]


def update_model(old_model: ModelPrev) -> ModelNew:
    """Translate a Maud old model to a new model."""
    reactions = []
    enzymes = []
    enzyme_reactions = []
    allosteries = []
    comp_inhibitions = []
    for reac in old_model.reaction:
        mechanism = (
            md.ReactionMechanism.REVERSIBLE_MICHAELIS_MENTEN
            if reac.reaction_mechanism.startswith("reversible")
            else md.ReactionMechanism.IRREVERSIBLE_MICHAELIS_MENTEN
        )
        reac_id = reac.id.replace("_", "")
        reactions.append(
            md.Reaction(
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
            enzymes.append(md.Enzyme(id=enz_id, name=enz.name, subunits=enz.subunits))
            enzyme_reactions.append(
                md.EnzymeReaction(enzyme_id=enz_id, reaction_id=reac_id)
            )
            if enz.modifier is not None:
                for modifier in enz.modifier:
                    if modifier.modifier_type == "competitive_inhibitor":
                        comp_inhibitions.append(
                            md.CompetitiveInhibition(
                                enzyme_id=enz_id,
                                reaction_id=reac_id,
                                metabolite_id=modifier.mic_id,
                                compartment_id="c",
                            )
                        )
                    else:
                        mod_type = (
                            md.ModificationType.INHIBITION
                            if modifier.modifier_type == "allosteric_inhibitor"
                            else md.ModificationType.ACTIVATION
                        )
                        allosteries.append(
                            md.Allostery(
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
            md.Metabolite(id=met_id, name=met.name, inchi_key=met.inchi_key)
        )
        comp_metabolites.append(
            md.MetaboliteInCompartment(
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
    """Serialize a model into toml.

    The intermediate JSON step is required so that dataclasses are
    serialized properly.
    """
    with open(out_file, "w") as f:
        toml.dump(json.loads(model.json()), f)


def read_old_maud(toml_file: str):
    with open(toml_file) as f:
        data = ModelPrev.parse_obj(toml.load(f))
    return data


@click.command()
@click.argument("old_toml", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path(dir_okay=False))
def cli_entry(old_toml: click.Path, output: click.Path):
    data = read_old_maud(old_toml)
    new_data = update_model(data)
    print(new_data.json())
    write_new_model(new_data, output)


if __name__ == "__main__":
    cli_entry()
