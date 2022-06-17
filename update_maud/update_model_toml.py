"""Script to transform Maud model to the current version."""
import json
import re
from typing import Optional

import click
import maud.data_model.kinetic_model as md
import toml
from pydantic import BaseModel, Field


UNDER_PAT = re.compile(r"_(?![a-z]$)")


class Compartment(BaseModel):
    id: str
    name: str
    volume: float


class MetabolitePrev(BaseModel):
    id: str = Field(alias="metabolite")
    compartment: str
    balanced: bool
    name: Optional[str] = None
    metabolite_inchi_key: Optional[str] = None


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
    mechanism: Optional[str] = "reversible_modular_rate_law"
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
    drain: Optional[list[DrainPrev]] = None

    class Config:
        allow_population_by_field_name = True


class ModelNew(BaseModel):
    compartment: list[Compartment]
    reaction: list[md.Reaction]
    enzyme: list[md.Enzyme]
    enzyme_reaction: list[md.EnzymeReaction]
    metabolite: list[md.Metabolite]
    metabolite_in_compartment: list[md.MetaboliteInCompartment] = Field(
        list[md.MetaboliteInCompartment], exclude={"__all__": {"id"}}
    )
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
            if reac.mechanism.startswith("reversible")
            else md.ReactionMechanism.IRREVERSIBLE_MICHAELIS_MENTEN
        )
        reac_id = reac.id.replace("_", "")
        reactions.append(
            md.Reaction(
                id=reac_id,
                name=reac.name,
                mechanism=mechanism,
                stoichiometry={
                    UNDER_PAT.sub("", k): v for k, v in reac.stoichiometry.items()
                },
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
                                metabolite_id=modifier.mic_id[:-2].replace("_", ""),
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
                                metabolite_id=modifier.mic_id[:-2].replace("_", ""),
                                compartment_id="c",
                                modification_type=mod_type,
                            )
                        )
    metabolites = []
    comp_metabolites = []
    for met in old_model.metabolite:
        met_id = met.id.replace("_", "")
        if not any(met_id == met.id for met in metabolites):
            metabolites.append(
                md.Metabolite(
                    id=met_id, name=met.name, inchi_key=met.metabolite_inchi_key
                )
            )
        comp_metabolites.append(
            md.MetaboliteInCompartment(
                metabolite_id=met_id,
                compartment_id=met.compartment,
                balanced=met.balanced,
            )
        )
    # TODO(jorge): not sure how phosphorylation looks lik
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


def exclude_id(model: BaseModel, keys: list[str]):
    for key in keys:
        for i in range(len(model[key])):
            del model[key][i]["id"]


def reaction_mech_to_name(model: dict):
    for reac in model["reaction"]:
        reac["mechanism"] = md.ReactionMechanism(reac["mechanism"]).name
    for allostery in model["allostery"]:
        allostery["modification_type"] = md.ModificationType(
            allostery["modification_type"]
        ).name


def write_new_model(model: ModelNew, out_file: str):
    """Serialize a model into toml.

    The intermediate JSON step is required so that dataclasses are
    serialized properly.
    """
    with open(out_file, "w") as f:
        model_dict = json.loads(model.json())
        exclude_id(
            model_dict,
            [
                "enzyme_reaction",
                "metabolite_in_compartment",
                "allostery",
                "competitive_inhibition",
            ],
        )
        reaction_mech_to_name(model_dict)
        toml.dump(model_dict, f)


def read_old_maud(toml_file: str):
    with open(toml_file) as f:
        data = ModelPrev.parse_obj(toml.load(f))
    return data


def update_model_toml(old_toml: click.Path, output: click.Path):
    data = read_old_maud(old_toml)
    new_data = update_model(data)
    write_new_model(new_data, output)
    return new_data
