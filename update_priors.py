"""Update priors file to newer maud impl."""

import sys

import click
import pandas as pd
import toml

from update_model_toml import UNDER_PAT, ModelNew

NEW_PRIOR_COLS = ["parameter", "row_id", "col_id", "location", "scale", "pct1", "pct99"]


def query_lookup(
    query_id: str,
    model: ModelNew,
    query_ns: str = "enzyme_id",
    lookup: str = "enzyme_reaction",
):
    """Lookup the corresponding reaction of an enzyme.

    Trivial since there were not promiscuous enzymes before.
    """
    enz_id = query_id.replace("_", "")
    return next(
        enz_reac.id
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
    df["col_id"] = None
    df = df.rename({"parameter_type": "parameter"}, axis=1)
    df.enzyme_id = df.enzyme_id.str.replace("_", "")
    df.mic_id = df.mic_id.apply(
        lambda x: UNDER_PAT.sub("", x) if isinstance(x, str) else x
    )
    df.experiment_id = df.experiment_id.str.replace("_", "")
    # make sure that every enzyme is actually in the final config
    df = df.loc[
        ~df.parameter.isin(["kcat", "conc_enzyme"])
        | df.enzyme_id.apply(lambda x: is_valid_enzyme(x, model)),
        :,
    ]
    # kcats
    df.loc[df.parameter == "kcat", "row_id"] = df.loc[
        df.parameter == "kcat", "enzyme_id"
    ].apply(lambda x: query_lookup(x, model))
    # kms
    kms = df[df.parameter == "km"]
    df.loc[df.parameter == "km", "row_id"] = kms.enzyme_id + "_" + kms.mic_id
    # ki
    df.loc[df.parameter == "ki", "row_id"] = (
        df.loc[df.parameter == "ki", "enzyme_id"].apply(
            lambda x: query_lookup(x, model)
        )
        + "_"
        + df.loc[df.parameter == "ki", "mic_id"]
    )
    # transfer and dissociation constants; phosphorilation
    dissociation_constants = df[df.parameter == "dissociation_constant"]
    df.loc[df.parameter == "dissociation_constant", "row_id"] = (
        dissociation_constants.enzyme_id + "_" + dissociation_constants.mic_id
    )
    df.loc[df.parameter.isin(["transfer_constant", "kcat_phos"]), "row_id"] = df.loc[
        df.parameter.isin(["transfer_constant", "kcat_phos"]), "enzyme_id"
    ]
    conc_phos = df[df.parameter == "conc_phos"]
    df.loc[df.parameter == "conc_phos", "row_id"] = (
        conc_phos.experiment_id + "_" + conc_phos.enzyme_id
    )
    # psi
    df.loc[df.parameter == "psi", "row_id"] = df.loc[
        df.parameter == "psi", "experiment_id"
    ]
    # conc_unbalanced
    df.loc[df.parameter == "conc_unbalanced", "row_id"] = df.loc[
        df.parameter == "conc_unbalanced", "experiment_id"
    ]
    df.loc[df.parameter == "conc_unbalanced", "col_id"] = df.loc[
        df.parameter == "conc_unbalanced", "mic_id"
    ]
    # protein concentrations
    df.loc[df.parameter == "conc_enzyme", "row_id"] = df.loc[
        df.parameter == "conc_enzyme", "experiment_id"
    ]
    df.loc[df.parameter == "conc_enzyme", "col_id"] = df.loc[
        df.parameter == "conc_enzyme", "enzyme_id"
    ]
    # drains
    drains = df[df.parameter == "drain"]
    df.loc[df.parameter == "drain", "row_id"] = (
        drains.experiment_id + "_" + drains.drain_id
    )
    df = df[NEW_PRIOR_COLS]
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
