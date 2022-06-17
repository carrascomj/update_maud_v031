"""Update priors file to newer maud impl."""

import re
from copy import deepcopy

import click
import pandas as pd
import toml

from .update_model_toml import UNDER_PAT, ModelNew

NEW_PRIOR_COLS = [
    "parameter",
    "metabolite",
    "compartment",
    "enzyme",
    "reaction",
    "experiment",
    "location",
    "scale",
    "pct1",
    "pct99",
]
# the comparment is represented by a single letter,
# prefixed by underscored at the end of the entity. 
COMP_PAT = re.compile(r"(.*)_([a-z]$)")


def query_lookup(
    query_id: str,
    model: ModelNew,
    query_ns: str = "enzyme_id",
    lookup: str = "enzyme_reaction",
    key: str = "reaction_id",
):
    """Lookup the corresponding reaction of an enzyme.

    Trivial since there were not promiscuous enzymes before.
    """
    enz_id = query_id.replace("_", "")
    return next(
        enz_reac.__getattribute__(key)
        for enz_reac in model.__getattribute__(lookup)
        if enz_reac.__getattribute__(query_ns) == enz_id
    )


def is_valid_enzyme(enzyme_id: str, model: ModelNew):
    if not isinstance(enzyme_id, str):
        return False
    enz_id = enzyme_id.replace("_", "")
    return any(
        enz_reac.reaction_id
        for enz_reac in model.enzyme_reaction
        if enz_reac.enzyme_id == enz_id
    )


def update_priors(priors_df: pd.DataFrame, model: ModelNew) -> pd.DataFrame:
    df = priors_df.copy()
    target_cols = deepcopy(NEW_PRIOR_COLS)
    df = df.rename(
        {
            "parameter_type": "parameter",
            "enzyme_id": "enzyme",
            "experiment_id": "experiment",
            "mic_id": "metabolite",
            # will fill out the other reactions later
            "drain_id": "reaction",
        },
        axis=1,
    )
    df["compartment"] = None
    # comparment separation of metabolites
    df.loc[~pd.isna(df.metabolite), "compartment"] = df.loc[
        ~pd.isna(df.metabolite), "metabolite"
    ].apply(lambda x: COMP_PAT.sub(r"\2", x))
    df.loc[~pd.isna(df.metabolite), "metabolite"] = df.loc[
        ~pd.isna(df.metabolite), "metabolite"
    ].apply(lambda x: COMP_PAT.sub(r"\1", x))
    # underscores are forbidden
    df.enzyme = df.enzyme.str.replace("_", "")
    df.experiment = df.experiment.str.replace("_", "")
    df.metabolite = df.metabolite.str.replace("_", "")
    df.parameter = df.parameter.str.replace("diss_t", "dissociation_constant")
    df.parameter = df.parameter.str.replace("conc_phos", "conc_pme")
    # make sure that every enzyme is actually in the final config
    df = df.loc[
        ~df.parameter.isin(["kcat", "conc_enzyme"])
        | df.enzyme.apply(lambda x: is_valid_enzyme(x, model)),
        :,
    ]
    # add reaction ids
    df.loc[df.parameter.isin(["kcat", "ki"]), "reaction"] = df.loc[
        df.parameter.isin(["kcat", "ki"]), "enzyme"
    ].apply(lambda x: query_lookup(x, model))
    if (df.parameter == "conc_phos").any():
        target_cols.append("phosphorylation_modifying_enzyme")
        # TODO(jorge): I have to look up an example of this
        raise NotImplementedError("Phosphorylation update is not implemented")
    # the previous steps are only succesful if this works
    df = df[target_cols]
    return df


@click.command()
@click.argument("old_priors", type=click.Path(exists=True, dir_okay=False))
@click.argument("new_toml", type=click.Path(exists=True, dir_okay=False))
@click.argument("output", type=click.Path(dir_okay=False))
def cli_entry(old_priors: click.Path, new_toml: click.Path, output: click.Path):
    with open(new_toml) as f:
        model = toml.load(f)
        model = ModelNew.parse_obj(model)
    update_priors(pd.read_csv(old_priors), model).to_csv(output, index=False)
