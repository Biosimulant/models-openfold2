# OpenFold2 validation evidence

BIO-006, 2026-09-20. Local contract verification passed; native inference and scientific qualification did not pass.

The individually authored adapter has 65 passing pytest checks on Python 3.10.20 with biosimulant 0.0.34. These cover manifest/runtime/Lab port parity, one-shot execution, request validation, ranked output preservation, residue and atom mapping, optional confidence fields, malformed outputs, bounded HTTP polling, failure handling, exact-request resume and credential absence. All prediction outputs in these tests are synthetic fixtures. No successful NVIDIA prediction has been captured for this entry.

Three authorized direct NVIDIA requests used the complete 46-residue 1CRN chain A query and the original 101-row MSA Search A3M. Two used identical sequence/MSA inputs; a distinct third case also supplied the complete public 1CRN mmCIF as an explicit template. This template case is a technical interface check and is not a held-out accuracy benchmark.

| Case | NVIDIA request ID | Request SHA-256 | Observed outcome |
|---|---|---|---|
| openfold2-native-reference-001 | 554155a3-741b-4dba-8a94-bea774f95b26 | 6ca3a3653a9a80a563feddb8e441705e1da7190e54609b0938b8502f871960a0 | HTTP 500; CUDA illegal memory access |
| openfold2-native-reference-002 | 269737f7-9998-42af-bc8b-daade1f32cb6 | 6ca3a3653a9a80a563feddb8e441705e1da7190e54609b0938b8502f871960a0 | HTTP 500; CUDA illegal memory access |
| openfold2-template-native-001 | af109cba-58c0-400d-b045-821f6cf88d7f | 37b349eb4690744765e3f2e8b9a7598ffa1081ef8dba4eca83488ef8fef81bfd | HTTP 500; CUDA illegal memory access |

All three returned the same 528-byte error response, SHA-256 `5793178b2e1443d1594fa2275ee3c7e0adca4fe102b106fda3c9abf9c34cf5e3`. They are terminal failed attempts; there are no active request handles. No structures, confidence tensors, scientific results or Hub run IDs were obtained. Identical requests are not being automatically retried.

The raw request and response bytes, receipts, fixture provenance, Lab input records and inspection reports are retained under the correspondingly named directories in the tracking output folder. An independent evidence audit verifies their hashes, MSA identity/columns and template coordinates. A separate request-parity check imports this adapter and rebuilds each request byte-for-byte from the recorded Lab inputs; all three match. Request parity is not independent verification of the prediction algorithm.

MSA artifact SHA-256: `3e92fc74f6daa19c0845e72f1b4ba3839d126162555852f58ec15ad74b2f404b` (10,637 bytes). Template artifact SHA-256: `23787562c427d7c1abe5420e86d5f1d0a6c7007dec1e8ce85645a6d69c32e8ba` (69,506 bytes), from https://files.rcsb.org/download/1CRN.cif. Adapter source SHA-256 for the audit: `18a488ba64da60a4fe894d2d723317453ad00585ed95216219af1674be92f42b`.

Remaining gates: provider service recovery or a verified supported exact-release endpoint; successful native and adapter inference; exact deployed code/weights/container identity and rights; justified structural benchmark criteria with leakage control; supported compatibility checkers/adapters; managed credential and runtime verification; independent inspection of captured prediction outputs; private review and exact public-release approval. No tolerance or scientific success is inferred from synthetic tests, one historical protein, model confidence or provider benchmark claims.
