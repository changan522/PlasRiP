"""Custom exception hierarchy for PlasRiP v2."""

from __future__ import annotations


class PlasRipError(Exception):
    """Base class for all PlasRiP-specific errors."""


class PlasRipConfigError(PlasRipError):
    """Raised when configuration or argument values are invalid."""


class PlasRipDependencyError(PlasRipError):
    """Raised when an external dependency is missing or unusable."""


class PlasRipInputError(PlasRipError):
    """Raised when the user input files are missing or malformed."""


class PlasRipModelError(PlasRipError):
    """Raised when model files or model loading fail."""


class PlasRipDataError(PlasRipError):
    """Raised when bundled or reference data are missing or malformed."""
