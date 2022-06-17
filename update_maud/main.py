"""Script to transform Maud model to the current version."""
import os
import shutil
from pathlib import Path

import click
import pandas as pd
import toml
from maud.data_model.maud_config import ODEConfig

from .update_model_toml import update_model_toml
from .update_priors import update_priors


def rename_keys(dict_: dict, key_map: dict) -> dict:
    """Update the keys in a dictionary given an {old_key: new_key} map."""
    return {
        key_map[old_key] if old_key in key_map else old_key: value
        for old_key, value in dict_.items()
    }


def update_config(old_config: Path, new_config: Path):
    with open(old_config) as f:
        config = toml.load(f)
    config = rename_keys(
        config,
        {
            "priors": "priors_file",
            "measurements": "measurements_file",
            "biological_config": "experimental_setup_file",
            "kinetic_model": "kinetic_model_file",
        },
    )
    config["ode_config"] = rename_keys(
        config["ode_config"],
        {
            "rel_tol_forward": "rel_tol",
            "abs_tol_forward": "abs_tol",
        },
    )
    allowed = ODEConfig.__dict__.keys()
    config["ode_config"] = {
        k: v for k, v in config["ode_config"].items() if k in allowed
    }
    with open(new_config, "w") as f:
        toml.dump(config, f)


def update_measurements(old_measurements: Path, new_measurements: Path):
    df = pd.read_csv(old_measurements)
    df.experiment_id = df.experiment_id.str.replace("_", "")
    df.loc[df.measurement_type == "flux", "target_id"] = df.loc[
        df.measurement_type == "flux", "target_id"
    ].str.replace("_", "")
    df.to_csv(new_measurements, index=False)


def update_biological_config(old_toml: Path, new_toml: Path):
    with open(old_toml) as f:
        config = toml.load(f)
    for exp in config["experiment"]:
        exp.update({"id": exp["id"].replace("_", "")})
    config["experiment"] = [
        rename_keys(exp, {"sample": "is_train", "predict": "is_test"})
        for exp in config["experiment"]
    ]
    with open(new_toml, "w") as f:
        toml.dump(config, f)


def update_dgf(config, data_path: Path, out_path: Path):
    """Remove all underscores in identifiers in dgf files."""
    mean_path = config["dgf_mean_file"]
    cov_path = config["dgf_covariance_file"]
    means = pd.read_csv(data_path / mean_path)
    cov = pd.read_csv(data_path / cov_path)
    means.metabolite = means.metabolite.str.replace("_", "")
    cov.metabolite = cov.metabolite.str.replace("_", "")
    cov.columns = cov.columns.str.replace("_", "")
    means.to_csv(out_path / mean_path, index=False)
    cov.to_csv(out_path / cov_path, index=False)


def copy_unaffected_files(config: dict[str, str], data_path: Path, out_path: Path):
    """Copy to new location the files that are not affected by the update."""
    if "user_inits_file" in config:
        shutil.copy(
            data_path / config["user_inits_file"], out_path / config["user_inits_file"]
        )


@click.command()
@click.argument("data_dir", type=click.Path(exists=True, dir_okay=True))
@click.argument("outdir", type=click.Path(dir_okay=True))
def cli_entry(data_dir: click.Path, outdir: click.Path):
    """Update an old maud data directory."""
    if not os.path.exists(outdir):
        os.mkdir(outdir)
    data_path = Path(data_dir)
    out_path = Path(outdir)
    with open(data_path / "config.toml") as f:
        config = toml.load(f)
    old_toml = data_path / config["kinetic_model"]
    new_toml = out_path / config["kinetic_model"]
    old_priors = data_path / config["priors"]
    new_priors = out_path / config["priors"]
    measurements_file = config["measurements"]
    bio_config_file = config["biological_config"]
    model = update_model_toml(old_toml, new_toml)
    update_priors(pd.read_csv(old_priors), model).to_csv(new_priors, index=False)
    update_measurements(data_path / measurements_file, out_path / measurements_file)
    update_biological_config(data_path / bio_config_file, out_path / bio_config_file)
    update_config(data_path / "config.toml", out_path / "config.toml")
    if "dgf_mean_file" in config:
        update_dgf(config, data_path, out_path)
    copy_unaffected_files(config, data_path, out_path)


if __name__ == "__main__":
    cli_entry()
