# PlasRiP v0.3.0

This release migrates the Mash-based cluster interface from PTU to MPC and reads the MPC metadata in Reference Database v3.0 (schema 2.0).

The new assignment outputs are `mpc.tsv` and `mpc.tsv.hits.tsv`. CLI options and Python API parameters use MPC names; `run_mpc_assignment` replaces `run_ptu_assignment`. Existing integrations must update field names and options.

Cluster identifiers retain their numeric values as `MPC000001` etc. The algorithm, Mash parameters, reference paths, model weights and application score-gap default (0.20) are unchanged. This interface update does not establish equivalence with the published plasmid taxonomic unit classification.

The update includes a source patch and installable wheel. The conda recipe is updated; an older conda binary is not a v0.3.0 build.
