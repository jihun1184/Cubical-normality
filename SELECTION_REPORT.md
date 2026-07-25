# Code-selection report

## Included as the public core

- The latest canonical `check_embedded_cubical_normality.py` implementation (552 lines), renamed to `src/cubical_normality/checker.py` without changing its mathematical logic.
- The latest canonical regression suite (445 lines; 88 tests in the inspected environment).
- The canonical Pipeline 1 v1 specification.
- A small packaging layer only: import surface, JSON CLI, examples, and Python project metadata.

## Why these files are the core

They directly implement the paper-facing decision procedure, are self-contained apart from the Python standard library, produce structured failure certificates, and have a deterministic regression suite. They are also the files explicitly retained in the session-41 canonical bundle.

## Excluded from the main repository

The full archive contained 708 Python file instances but only 244 distinct code blobs. Most excluded files fall into these categories:

- repeated copies carried across handoff ZIPs;
- exploratory witness searches and large random sweeps;
- session-specific counterfactual, signature-9, reachability, and Stage 3 analyses;
- scripts coupled to local JSON/CSV/checkpoint files or undocumented directory layouts;
- superseded checker/test versions;
- one-off diagnostics whose conclusions were absorbed into the proof or canonical specification.

These files are valuable provenance, but publishing them beside the core checker would make the repository harder to audit and would blur the boundary between the theorem-facing implementation and historical exploratory computation.

## Optional later release

A separate `research-archive` branch or Zenodo supplement could preserve the exploratory scripts together with their exact input data, environment lockfile, command manifest, and expected output hashes. They should not be mixed into the minimal software artifact until those dependencies are normalized.
