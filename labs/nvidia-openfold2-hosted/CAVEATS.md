# OpenFold2 caveats

Current native prediction is blocked by NVIDIA CUDA illegal-memory-access errors on the complete public 1CRN query with its real MSA. Repeating the identical request returned the same HTTP 500 under a new provider ID. Neither is a successful inference run. Local code/tests cannot repair the hosted GPU process or establish its hidden checkpoint identity.

Current NVIDIA model card/NIM target is 2.6.0; hosted published API remains 2.1.0 and contains legacy fields. HHR templates and relaxation are unsupported in current NIM. Model-card 1000-residue and current hardware 2048-residue bounds differ; the hosted adapter uses the published 1000-residue boundary. Exact deployed image/weights, managed runtime and secure credentials remain unverified.

MSA data comes from a genuine NVIDIA MSA Search request with explicit database name and captured bytes, but the database/index digest is unavailable. 1CRN is an old public structure and may be present in training or MSA databases. Supplying the same 1CRN structure as a template tests the input interface, not independent prediction accuracy.

Confidence, pLDDT, PAE and pTM are model estimates with different scales and identities. Confidence is not observed accuracy. A PDB B-factor is not automatically an experimental uncertainty estimate or a separately verified pLDDT output. Proposed compatibility profiles are unreleased; scientific thresholds and independently reviewed benchmarks remain blockers.

No Hub run, private release, public qualification or clinical-use claim has been established. No model weights, container or runtime were downloaded/installed.
