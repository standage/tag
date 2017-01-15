# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
### Added
- An index class for efficient in-memory access of sequence features.
- Module for mRNA handling, with a function for selecting the primary mRNA from
  a gene or other feature.
- New CLI command `tag pmrna`.
- A new `Score` class for internal handling of feature scores. Not yet included
  in the API, and may not ever be.

### Changed
- Modules focused on classes / data structure now support more concise imports
  (for example, `from tag import Feature` and `tag.Feature` now supported and
  preferred over `from tag.feature import feature` and `tag.feature.Feature`).

### Fixed
- Resolved a bug with the GFF3Writer failing to print `##FASTA` directives
  before writing sequences to output.

## [0.1.1] - 2016-12-19
### Changed
- CLI implemented using `entry_points` instead of a dedicated script.

### Fixed
- Entry type inference now correct by inheriting from `object`.

## [0.1.0] - 2016-12-16
### Added
- Basic data structures
    - Range
    - Comment
    - Directive
    - Sequence
    - Feature
- Annotation I/O
    - GFF3Reader
    - GFF3Writer
- Composable generator functions for streaming annotation processing
- A command line interface through the `tag` script
- Package scaffolding
    - README
    - documentation
    - license
    - changelog
    - various config files
