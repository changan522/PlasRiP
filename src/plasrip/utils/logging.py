"""Simple logger factory for PlasRiP."""

from __future__ import annotations

import logging
from pathlib import Path


def setup_logger(name: str, log_file: str | None = None, quiet: bool = False) -> logging.Logger:
    """Create or reconfigure a logger with stream and optional file handlers."""
    logger = logging.getLogger(name)
    logger.setLevel(logging.INFO)
    logger.handlers.clear()
    logger.propagate = False

    formatter = logging.Formatter("[%(levelname)s] %(message)s")

    if not quiet:
        sh = logging.StreamHandler()
        sh.setFormatter(formatter)
        logger.addHandler(sh)

    if log_file is not None:
        path = Path(log_file).resolve()
        path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(path, encoding="utf-8")
        fh.setFormatter(formatter)
        logger.addHandler(fh)

    return logger


def get_logger(name: str = "plasrip") -> logging.Logger:
    """Return a logger by name."""
    return logging.getLogger(name)


def close_logger(logger: logging.Logger) -> None:
    """Close all handlers attached to a logger."""
    handlers = list(logger.handlers)
    for handler in handlers:
        handler.close()
        logger.removeHandler(handler)
