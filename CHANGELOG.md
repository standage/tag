# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).

## [Unreleased]
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
