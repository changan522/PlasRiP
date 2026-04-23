"""Console helpers for old-style colorful per-plasmid progress output."""

from __future__ import annotations

import os
import sys
from pathlib import Path

RESET = "\033[0m"
BOLD = "\033[1m"
CYAN = "\033[96m"
GREEN = "\033[92m"
YELLOW = "\033[93m"
RED = "\033[91m"
BLUE = "\033[94m"

HEADER_ART = (
    "\n"
    "     *            *               *                   *\n"
    " *        *                *             *       *     *\n"
    f"{CYAN}||{RESET}{BOLD}{RED} ____  _           ____  _ ____  {RESET}{CYAN}||{RESET}\n"
    f"{CYAN}||{RESET}{BOLD}{RED}|  _ \\| | __ _ ___|  _ \\(_)  _ \\ {RESET}{CYAN}||{RESET}\n"
    f"{CYAN}||{RESET}{BOLD}{RED}| |_) | |/ _` / __| |_) | | |_) |{RESET}{CYAN}||{RESET}\n"
    f"{CYAN}||{RESET}{BOLD}{RED}|  __/| | (_| \\__ \\  _ <| |  __/ {RESET}{CYAN}||{RESET}\n"
    f"{CYAN}||{RESET}{BOLD}{RED}|_|   |_|\\__,_|___/_| \\_\\_|_|    {RESET}{CYAN}||{RESET}\n"
    " *        *           *        *      *          *       \n"
    "     *         *               *              *          \n"
    "\n"
)


def _write(msg: str) -> None:
    sys.stderr.write(msg + "\n")
    sys.stderr.flush()


def print_header() -> None:
    _write(HEADER_ART.rstrip("\n"))


def banner(title: str) -> None:
    line = "=" * 70
    _write(line)
    _write(f"{title:^{70}}")
    _write(line)


def print_step(plasmid: str, step_name: str) -> None:
    _write(f"{BLUE}[STEP]{RESET} {plasmid} :: {step_name}")


def print_info(message: str) -> None:
    _write(f"{CYAN}[INFO]{RESET} {message}")


def print_done(message: str) -> None:
    _write(f"{GREEN}[DONE]{RESET} {message}")


def print_warn(message: str) -> None:
    _write(f"{YELLOW}[WARN]{RESET} {message}")


def print_error(message: str) -> None:
    _write(f"{RED}[ERROR]{RESET} {message}")


def print_plasmid_header(index: int, total: int, fasta_path: str) -> None:
    name = os.path.basename(str(fasta_path))
    line = "-" * 70
    _write(line)
    _write(f"[{index}/{total}] {name}")
    _write(line)


def print_final_summary(total: int, completed: int, failed: int, out_dir: str | Path | None = None) -> None:
    line = "=" * 40
    _write(line)
    if out_dir is not None:
        _write(f"[DONE] output dir       : {Path(out_dir).resolve()}")
    _write(f"[DONE] total plasmids   : {total}")
    _write(f"[DONE] completed        : {completed}")
    _write(f"[DONE] failed           : {failed}")
    _write(line)
