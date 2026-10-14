"""Public Python API for PlasRiP v0.3.

The analysis functions accept a database directory (or ``PLASRIP_DATABASE``).
For applications that need to inspect or validate the external reference
database directly, the same schema-2.0 database interface used by the CLI is
exposed here as a small public API.
"""

from __future__ import annotations

from os import PathLike

from .annotate import run_annotation
from .assign import run_mpc_assignment
from .core import run_core_analysis
from .database import (
    PlasRipReferenceDatabase,
    load_database,
    validate_database_dir,
)
from .pipeline import run_pipeline
from .predict import run_host_prediction


def load_reference_database(
    database_dir: str | PathLike[str] | None = None,
    *,
    check_sketches: bool = True,
) -> PlasRipReferenceDatabase:
    """Validate and load a PlasRiP schema-2.0 reference database.

    ``database_dir`` may be omitted when ``PLASRIP_DATABASE`` is set.
    """
    return load_database(database_dir, check_sketches=check_sketches)


def validate_reference_database(
    database_dir: str | PathLike[str] | None = None,
    *,
    check_sketches: bool = False,
    check_all_sketches: bool = False,
) -> dict[str, object]:
    """Return the same database-validation report used by the CLI.

    Set ``check_sketches`` to validate MPC representative sketches, or
    ``check_all_sketches`` for full release QC of the Mash index.
    """
    return validate_database_dir(
        database_dir,
        check_sketches=check_sketches,
        check_all_sketches=check_all_sketches,
    )


__all__ = [
    "PlasRipReferenceDatabase",
    "load_reference_database",
    "validate_reference_database",
    "run_host_prediction",
    "run_mpc_assignment",
    "run_core_analysis",
    "run_annotation",
    "run_pipeline",
]
