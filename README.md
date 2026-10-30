# ActivityFinder analysis workspace

This workspace contains two separate analysis projects. Their code, data requirements,
and setup are independent; use the README in the relevant project directory for details.

## Projects

- **[Comparison_analysis](Comparison_analysis/README.md)** — compares ActivityFinder
  results with BioChemGraph and contains the related sequence, molecule, mutation,
  and visualization analyses.
- **[ActivityDB_analysis](ActivityDB_analysis/README.md)** — a separate ActivityDB
  analysis project. Its files are maintained independently of `Comparison_analysis`.

## Working in this workspace

Run commands from the project directory they belong to, and keep project-specific
inputs, outputs, and environment configuration with that project. In particular,
the Comparison project expects external data and database access that are not
provided by this workspace; see its README before running its scripts.

The two directories are separate codebases even though they are checked out under
the same workspace. Changes to one project should not be assumed to apply to the
other.
