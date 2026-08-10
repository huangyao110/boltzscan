"""Adapter between vendored upstream Boltz and BoltzScan ODE variants.

The wheel carries an unmodified, pinned Boltz inference source tree. Standard
Boltz1/2 therefore use their native CLI defaults. This worker implements the
project-owned ODE protocols without forking that CLI: it maps the public ODE
name to its base family and replaces only that family's diffusion-parameter
dataclass before entering the native Click command.
"""

import sys
from dataclasses import dataclass


@dataclass
class BoltzODEDiffusionParams:
    """BoltzScan probability-flow-style ODE sampling parameters."""

    gamma_0: float = 0.0
    gamma_min: float = 1.0
    noise_scale: float = 1.0
    rho: float = 7
    step_scale: float = 1.0
    sigma_min: float = 0.0001
    sigma_max: float = 160.0
    sigma_data: float = 16.0
    P_mean: float = -1.2
    P_std: float = 1.5
    coordinate_augmentation: bool = True
    alignment_reverse_diff: bool = True
    synchronize_sigmas: bool = True


def _configure_model(argv):
    """Patch the selected base family and return native Boltz CLI arguments."""
    from boltz import main as boltz_main

    argv = list(argv)
    try:
        model_index = argv.index("--model") + 1
        public_model = argv[model_index]
    except (ValueError, IndexError):
        return boltz_main, argv

    if public_model == "boltz1_ode":
        boltz_main.BoltzDiffusionParams = BoltzODEDiffusionParams
        argv[model_index] = "boltz1"
    elif public_model == "boltz2_ode":
        boltz_main.Boltz2DiffusionParams = BoltzODEDiffusionParams
        argv[model_index] = "boltz2"
    return boltz_main, argv


def main(argv=None):
    boltz_main, native_argv = _configure_model(
        sys.argv[1:] if argv is None else argv
    )
    return boltz_main.cli.main(
        args=native_argv,
        prog_name="boltz",
        standalone_mode=True,
    )


if __name__ == "__main__":
    main()
