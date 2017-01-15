# tag: Toolkit for Annotating Genomes

![Supported Python versions](https://img.shields.io/pypi/pyversions/tag.svg)
[![PyPI version](https://img.shields.io/pypi/v/tag.svg)](https://pypi.python.org/pypi/tag)
[![GenHub build status](https://img.shields.io/travis/standage/tag.svg)](https://travis-ci.org/standage/tag)
[![codecov.io coverage](https://img.shields.io/codecov/c/github/standage/tag.svg)](https://codecov.io/github/standage/tag)

> *Computational biology is 90% text formatting and ID cross-referencing!*  
> -- discouraged graduate students everywhere

**tag** is a free open-source software package for analyzing genome annotation data.

```python
# Compute number of exons per gene
import tag
reader = tag.GFF3Reader(infilename='/data/genomes/mybug.gff3.gz')
for gene in tag.select.features(reader, type='gene'):
    exons = [feat for feat in gene if feat.type == exon]
    print('num exons:', len(exons))
```

To install the most recent stable release execute `pip install tag` from your terminal.
Full installation instructions and project documentation are available at https://tag.readthedocs.io.
