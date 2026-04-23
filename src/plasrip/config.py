"""Central configuration defaults for PlasRiP v2 (3.2 unified)."""

from __future__ import annotations

from pathlib import Path

DEFAULT_THREADS = 1
DEFAULT_TOPK = 10
DEFAULT_PTU_THRESHOLD = 0.05
DEFAULT_AMR_DB = "ncbi"
DEFAULT_TAX_LEVEL = "species"

# 3.2 host-range semantics
DEFAULT_SCORE_DELTA = 0.20
DEFAULT_RANGE_DELTA_QUANTILE = 0.95

# 3.2 model packaging
DEFAULT_MODEL_FILENAME = "two_tower_best_3.2.pth"
REQUIRED_MODEL_FILENAMES = ["two_tower_best_3.2.pth"]

DEFAULT_LOG_NAME = "PlasRip.log"
DEFAULT_PREDICT_LOG_NAME = "PlasRip.predict.log"
DEFAULT_CORE_LOG_NAME = "PlasRip.core.log"
DEFAULT_ANNOTATE_LOG_NAME = "PlasRip.annotate.log"
DEFAULT_ASSIGN_LOG_NAME = "PlasRip.assign.log"
DEFAULT_REPLICON_LOG_NAME = "PlasRip.replicon.log"

DEFAULT_HOSTRANGE_OUT = "hostrange.tsv"
DEFAULT_HOSTRANGE_TOP_HITS_OUT = "hostrange.top_hits.tsv"
DEFAULT_HOSTRANGE_TAXONOMY_OUT = "hostrange.taxonomy.tsv"
DEFAULT_HOSTRANGE_SET_OUT = "hostrange.set.tsv"
DEFAULT_HOSTRANGE_ALL_SCORES_OUT = "hostrange.all_scores.tsv"

DEFAULT_PTU_OUT = "ptu.tsv"
DEFAULT_PTU_HITS_SUFFIX = ".hits.tsv"

DEFAULT_AMR_OUT = "amr.tsv"
DEFAULT_REPLICON_OUT = "plasmidfinder.tsv"
DEFAULT_MOBILITY_OUT = "mobility_summary.tsv"

DEFAULT_CORE_SUMMARY_OUT = "core_summary.tsv"
DEFAULT_ANNOTATION_SUMMARY_OUT = "annotation_summary.tsv"
DEFAULT_SUMMARY_OUT = "summary.tsv"

PTU_METADATA_FILENAME = "PTU-metadata.csv"
PTU_MEMBER_FILENAME = "PTUmember.csv"
PTU_INFO_FILENAME = "PTU_info.tsv"

HOSTRANGE_SUBDIR = "hostrange"
CONJSCAN_SUBDIR = "CONJScan"


def get_package_dir() -> Path:
    return Path(__file__).resolve().parent


def get_data_dir() -> Path:
    return get_package_dir() / "data"


def get_hostrange_data_dir() -> Path:
    return get_data_dir() / HOSTRANGE_SUBDIR


def get_conjscan_dir() -> Path:
    return get_data_dir() / CONJSCAN_SUBDIR


def get_ptu_metadata_path() -> Path:
    return get_data_dir() / PTU_METADATA_FILENAME


def get_ptu_member_path() -> Path:
    return get_data_dir() / PTU_MEMBER_FILENAME


def get_ptu_info_path() -> Path:
    return get_data_dir() / PTU_INFO_FILENAME
