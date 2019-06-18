# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased
### Added
- New `tag.select.merge` function to efficiently merge sorted feature streams (see #68).
- New `tag.locus.loci` function to efficiently determine locus coordinates from multiple annotation streams (see #71).

### Changed
- Minor updates to compensate for a couple years' worth of neglect (see #64).
- Refactored `GFF3Reader` to better support processing of unsorted GFF3 data (see #65).
- Modified `GFF3Writer` to intermittently output separator directives (`###`) in large stretches of simple (childless) features (see #68).
- Reorganized test suite code and data (see #70).



## [0.3.3] - 2017-06-21
### Fixed
- Missing `extent` query from the index implementation.
- Aliased `index.keys()` to `index.seqids`.

### Added
- Script `tag sum` to provide very basic summaries of genomic GFF3 files.

## [0.3.2] - 2017-04-19
### Added
- Pseudo-features for better handling and sorting of top-level multi-features.
- A new `primary_transcript` filter as a generalization of the `primary_mrna`
  function.
- A new function to query features for NCBI `GeneID` values.
- A new function to traverse a feature and all of its children to collect all
  attribute values associate with a given key.

### Fixed
- Bug with non-protein coding genes and the `tag.mrna.primary` filter (now
  `primary_mrna` in the `tag.transcript` module).
- Bug with how the GFF3 writer handles multi-feature IDs.

### Changed
- Refactored the `mrna` module, extended it, and renamed it to `transcript` to
  reflect its new and broader scope.

## [0.3.1] - 2017-02-02
### Fixed
- Range overlap queries accidentally left out of the previous release.

## [0.3.0] - 2017-02-02
### Added
- New convenience functions in the `Range` and `Feature` classes for range and
  point overlap queries.

## [0.2.0] - 2017-01-18
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
