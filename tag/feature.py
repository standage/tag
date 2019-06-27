#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of tag (http://github.com/standage/tag) and is licensed
# under the BSD 3-clause license: see LICENSE.
# -----------------------------------------------------------------------------

import tag
from tag.comment import Comment
from tag.directive import Directive
from tag.range import Range
from tag.sequence import Sequence
from tag.score import Score


class Feature(object):
    """Represents a feature entry from a GFF3 file.

    >>> gene = Feature('contig1', 'gene', 999, 7500, source='snap', strand='+',
    ...                attrstr='ID=gene1')
    >>> gene.seqid
    'contig1'
    >>> gene.source
    'snap'
    >>> gene.type
    'gene'
    >>> gene.range
    [999, 7500)
    >>> gene.score is None
    True
    >>> gene.strand
    '+'
    >>> gene.phase is None
    True
    >>> gene.get_attribute('ID')
    'gene1'
    >>> gene.get_attribute('NotARegisteredAttribute')
    >>> gene.slug
    'gene@contig1[1000, 7500]'
    """

    @staticmethod
    def from_gff3(data):
        fields = data.split('\t')
        assert len(fields) == 9
        if fields[6] not in ['+', '-', '.']:
            raise ValueError('invalid strand "{}"'.format(fields[6]))
        if fields[7] not in ['0', '1', '2', '.']:
            raise ValueError('invalid phase "{}"'.format(fields[7]))
        phase = None if fields[7] == '.' else int(fields[7])

        feat = Feature(
            fields[0], fields[2], int(fields[3]) - 1, int(fields[4]),
            source=fields[1], score=Score.from_str(fields[5]),
            strand=fields[6], phase=phase, attrstr=fields[8]
        )
        return feat

    def __init__(self, seqid, ftype, start, end, source='tag', score=None,
                 strand=None, phase=None, attrstr=None):
        # Core data
        self._seqid = seqid
        self._source = source
        self._type = ftype
        self._range = Range(start, end)
        self._score = score if isinstance(score, Score) else Score(score)
        self._strand = strand
        if strand not in ['-', '+', '.', None]:
            raise ValueError('invalid strand "{}"'.format(strand))
        self._phase = phase
        if phase not in [0, 1, 2, None]:
            raise ValueError('invalid phase "{}"'.format(phase))
        self._attrs = dict()
        if attrstr:
            self._attrs = self.parse_attributes(attrstr)

        # Ancillary data
        self.children = None
        self.multi_rep = None
        self.siblings = None
        self._pseudo = False

    def __str__(self):
        """String representation of the feature, sans children."""
        if self.is_pseudo:
            return ''

        phase = '.'
        if self.phase is not None:
            phase = str(self.phase)
        return '\t'.join([
            self.seqid, self.source, self.type, str(self.start + 1),
            str(self.end), str(self._score), self.strand, phase,
            self.attributes
        ])

    def __repr__(self):
        """Full representation of the feature, with children"""
        string = ''
        for feature in self:
            if string != '':
                string += '\n'
            string += str(feature)
        return string

    def __len__(self):
        return len(self._range)

    def like(self, other):
        if not isinstance(other, Feature):
            return False
        elif self.seqid != other.seqid or self._range != other._range:
            return False
        elif self.type != other.type or self.source != other.source:
            return False
        else:
            return True

    def __lt__(self, other):
        if isinstance(other, Directive) or isinstance(other, Comment):
            return False
        elif isinstance(other, Sequence):
            return True
        assert isinstance(other, Feature)
        if self.seqid != other.seqid:
            return self.seqid < other.seqid
        elif self._range != other._range:
            return self._range < other._range
        elif self.type != other.type:
            return self.type > other.type
        return self.source < other.source

    def __le__(self, other):
        if isinstance(other, Directive) or isinstance(other, Comment):
            return False
        elif isinstance(other, Sequence):
            return True
        assert isinstance(other, Feature)
        if self.seqid != other.seqid:
            return self.seqid < other.seqid
        elif self._range != other._range:
            return self._range < other._range
        elif self.type != other.type:
            return self.type > other.type
        return self.source <= other.source

    def __gt__(self, other):
        return not self.__le__(other)

    def __ge__(self, other):
        return not self.__lt__(other)

    def __iter__(self):
        """Generator iterates through a feature and all its subfeatures."""
        sorted_features = list()
        root = self
        if self.is_pseudo:
            root = self.children[0].multi_rep
        root._visit(L=sorted_features, marked={}, tempmarked={})
        for feat in sorted_features:
            yield feat

    def _visit(self, L, marked, tempmarked):
        """
        Sort features topologically.

        This recursive function uses depth-first search to find an ordering of
        the features in the feature graph that is sorted both topologically and
        with respect to genome coordinates.

        Implementation based on Wikipedia's description of the algorithm in
        Cormen's *Introduction to Algorithms*.
        http://en.wikipedia.org/wiki/Topological_sorting#Algorithms

        There are potentially many valid topological sorts of a feature graph,
        but only one that is also sorted with respect to genome coordinates
        (excluding different orderings of, for example, exons and CDS features
        with the same coordinates). Iterating through feature children in
        reversed order (in this functions' inner-most loop) seems to be the key
        to sorting with respect to genome coordinates.
        """
        assert not self.is_pseudo
        if self in tempmarked:
            raise ValueError('feature graph is cyclic')
        if self not in marked:
            tempmarked[self] = True
            features = list()
            if self.siblings is not None and self.is_toplevel:
                features.extend(reversed(self.siblings))
            if self.children is not None:
                features.extend(reversed(self.children))
            if len(features) > 0:
                for feature in features:
                    feature._visit(L, marked, tempmarked)
            marked[self] = True
            del tempmarked[self]
            L.insert(0, self)

    def add_child(self, child, rangecheck=False):
        """Add a child feature to this feature."""
        assert self.seqid == child.seqid, \
            (
                'seqid mismatch for feature {} ({} vs {})'.format(
                    self.fid, self.seqid, child.seqid
                )
            )
        if rangecheck is True:
            assert self._strand == child._strand, \
                ('child of feature {} has a different strand'.format(self.fid))
            assert self._range.contains(child._range), \
                (
                    'child of feature {} is not contained within its span '
                    '({}-{})'.format(self.fid, child.start, child.end)
                )
        if self.children is None:
            self.children = list()
        self.children.append(child)
        self.children.sort()

    @property
    def is_pseudo(self):
        return self._pseudo is True

    def pseudoify(self):
        """
        Derive a pseudo-feature parent from the given multi-feature.

        The provided multi-feature does not need to be the representative. The
        newly created pseudo-feature has the same seqid as the provided multi-
        feature, and spans its entire range. Otherwise, the pseudo-feature is
        empty. It is used only for convenience in sorting.
        """
        assert self.is_toplevel
        assert self.is_multi
        assert len(self.multi_rep.siblings) > 0
        rep = self.multi_rep

        start = min([s.start for s in rep.siblings + [rep]])
        end = max([s.end for s in rep.siblings + [rep]])

        parent = Feature(self._seqid, None, start, end, strand=self._strand)
        parent._pseudo = True
        for sibling in rep.siblings + [rep]:
            parent.add_child(sibling, rangecheck=True)
        parent.children = sorted(parent.children)
        rep.siblings = sorted(rep.siblings)

        return parent

    @property
    def num_children(self):
        if self.children is None:
            return 0
        return len(self.children)

    @property
    def fid(self):
        return self.get_attribute('ID')

    @property
    def slug(self):
        """
        A concise slug for this feature.

        Unlike the internal representation, which is 0-based half-open, the
        slug is a 1-based closed interval (a la GFF3).
        """
        return '{:s}@{:s}[{:d}, {:d}]'.format(self.type, self.seqid,
                                              self.start + 1, self.end)

    @property
    def is_multi(self):
        return self.multi_rep is not None

    @property
    def is_complex(self):
        return self.children is not None or self.is_multi

    @property
    def is_toplevel(self):
        if self.is_pseudo:
            return True
        return self.get_attribute('Parent') is None

    def add_sibling(self, sibling):
        """
        Designate this a multi-feature representative and add a co-feature.

        Some features exist discontinuously on the sequence, and therefore
        cannot be declared with a single GFF3 entry (which can encode only a
        single interval). The canonical encoding for these types of features is
        called a multi-feature, in which a single feature is declared on
        multiple lines with multiple entries all sharing the same feature type
        and ID attribute. This is commonly done with coding sequence (CDS)
        features.

        In this package, each multi-feature has a single "representative"
        feature object, and all other objects/entries associated with that
        multi-feature are attached to it as "siblings".

        Invoking this method will designate the calling feature as the
        multi-feature representative and add the argument as a sibling.
        """
        assert self.is_pseudo is False
        if self.siblings is None:
            self.siblings = list()
            self.multi_rep = self
        sibling.multi_rep = self
        self.siblings.append(sibling)

    @property
    def seqid(self):
        return self._seqid

    @seqid.setter
    def seqid(self, newseqid):
        """When modifying seqid, make sure to update seqid of children also."""
        for feature in self:
            feature._seqid = newseqid

    @property
    def source(self):
        return self._source

    @source.setter
    def source(self, newsource):
        """When modifying source, also update children with matching source."""
        oldsource = self.source
        for feature in self:
            if feature.source == oldsource:
                feature._source = newsource

    @property
    def type(self):
        if self.is_pseudo:
            return self.children[0].type
        return self._type

    @type.setter
    def type(self, newtype):
        """If the feature is a multifeature, update all entries."""
        self._type = newtype
        if self.is_multi:
            for sibling in self.multi_rep.siblings:
                sibling._type = newtype

    @property
    def start(self):
        return self._range.start

    @property
    def end(self):
        return self._range.end

    @property
    def range(self):
        return Range(self._range.start, self._range.end)

    def set_coord(self, start, end):
        """Manually reset the feature's coordinates."""
        self._range = Range(start, end)

    def transform(self, offset, newseqid=None):
        """Transform the feature's coordinates by the given offset."""
        for feature in self:
            feature._range.transform(offset)
            if newseqid is not None:
                feature.seqid = newseqid

    @property
    def score(self):
        return self._score.value

    @score.setter
    def score(self, newscore):
        self._score = Score(newscore)

    @property
    def strand(self):
        return '.' if self._strand is None else self._strand

    @strand.setter
    def strand(self, newstrand):
        if newstrand not in ['-', '+', '.', None]:
            raise ValueError('invalid strand "{}"'.format(newstrand))
        self._strand = newstrand if newstrand in ['-', '+'] else None

    @property
    def phase(self):
        return self._phase

    @property
    def attributes(self):
        if len(self._attrs) == 0:
            return '.'

        attrs = list()
        if 'ID' in self._attrs:
            attrs.append('ID=' + self.get_attribute('ID'))
        if 'Parent' in self._attrs:
            parent = self.get_attribute('Parent', as_string=True)
            attrs.append('Parent=' + parent)
        if 'Name' in self._attrs:
            name = self.get_attribute('Name', as_string=True)
            attrs.append('Name=' + name)
        for attrkey in sorted(self._attrs):
            if attrkey in ['ID', 'Parent', 'Name']:
                continue
            value = self.get_attribute(attrkey, as_string=True)
            if isinstance(value, float):
                attrs.append('{}={:.4f}'.format(attrkey, value))
            else:
                attrs.append('{}={}'.format(attrkey, value))
        return ';'.join(attrs)

    def add_attribute(self, attrkey, attrvalue, append=False, oldvalue=None):
        """
        Add an attribute to this feature.

        Feature attributes are stored as nested dictionaries.

        Each feature can only have one ID, so ID attribute mapping is 'string'
        to 'string'. All other attributes can have multiple values, so mapping
        is 'string' to 'dict of strings'.

        By default, adding an attribute that already exists will cause the old
        value to be overwritten. If the `append` option is true, the new
        attribute value will not overwrite the old value, but will be appended
        as a second value. (Note: ID attributes can have only 1 value.)

        If the `oldvalue` option is set, the new value will replace the old
        value. This is necessary for updating an attribute that has multiple
        values without completely overwriting all old values. (Note: The
        `append` option is ignored when `oldvalue` is set.)
        """
        # Handle ID/Parent relationships
        if attrkey == 'ID':
            if self.children is not None:
                oldid = self.get_attribute('ID')
                for child in self.children:
                    child.add_attribute('Parent', attrvalue,
                                        oldvalue=oldid)
            self._attrs[attrkey] = attrvalue
            if self.is_multi:
                self.multi_rep._attrs[attrkey] = attrvalue
                for sibling in self.multi_rep.siblings:
                    sibling._attrs[attrkey] = attrvalue
            return

        # Handle all other attribute types
        if oldvalue is not None:
            if attrkey in self._attrs:
                assert oldvalue in self._attrs[attrkey]
                del self._attrs[attrkey][oldvalue]
        if attrkey not in self._attrs or append is False:
            self._attrs[attrkey] = dict()
        self._attrs[attrkey][attrvalue] = True

    def get_attribute(self, attrkey, as_string=False, as_list=False):
        """
        Get the value of an attribute.

        By default, returns a string for ID and attributes with a single value,
        and a list of strings for attributes with multiple values. The
        `as_string` and `as_list` options can be used to force the function to
        return values as a string (comma-separated in case of multiple values)
        or a list.
        """
        assert not as_string or not as_list
        if attrkey not in self._attrs:
            return None
        if attrkey == 'ID':
            return self._attrs[attrkey]
        attrvalues = list(self._attrs[attrkey])
        attrvalues.sort()
        if len(attrvalues) == 1 and not as_list:
            return attrvalues[0]
        elif as_string:
            return ','.join(attrvalues)
        return attrvalues

    def drop_attribute(self, attrkey):
        """Drop the specified attribute from the feature."""
        if attrkey in self._attrs:
            del self._attrs[attrkey]

    def get_attribute_keys(self):
        """Return a list of all this feature's attribute keys."""
        return sorted(list(self._attrs))

    def parse_attributes(self, attrstring):
        """
        Parse an attribute string.

        Given a string with semicolon-separated key-value pairs, populate a
        dictionary with the given attributes.
        """
        if attrstring in [None, '', '.']:
            return dict()

        attributes = dict()
        keyvaluepairs = attrstring.split(';')
        for kvp in keyvaluepairs:
            if kvp == '':
                continue
            key, value = kvp.split('=')
            if key == 'ID':
                assert ',' not in value
                attributes[key] = value
                continue
            values = value.split(',')
            valdict = dict((val, True) for val in values)
            attributes[key] = valdict
        return attributes

    def attribute_crawl(self, key):
        """
        Grab all attribute values associated with the given feature.

        Traverse the given feature (and all of its descendants) to find all
        values associated with the given attribute key.

        >>> reader = tag.GFF3Reader(
        ...     tag.tests.data_stream('otau-no-seqreg.gff3')
        ... )
        >>> features = tag.select.features(reader)
        >>> for feature in features:
        ...     names = feature.attribute_crawl('Name')
        ...     print(sorted(list(names)))
        ['Ot01g00060', 'XM_003074019.1', 'XP_003074065.1']
        ['Ot01g00070', 'XM_003074020.1', 'XP_003074066.1']
        ['Ot01g00080', 'XM_003074021.1', 'XP_003074067.1']
        ['Ot01g00090', 'XM_003074022.1', 'XP_003074068.1']
        ['Ot01g00100', 'XM_003074023.1', 'XP_003074069.1']
        ['Ot01g00110', 'XM_003074024.1', 'XP_003074070.1']
        """
        union = set()
        for feature in self:
            values = feature.get_attribute(key, as_list=True)
            if values is not None:
                union.update(set(values))
        return union

    @property
    def ncbi_geneid(self):
        """
        Retrieve this feature's NCBI GeneID if it's present.

        NCBI GFF3 files contain gene IDs encoded in **Dbxref** attributes
        (example: `Dbxref=GeneID:103504972`). This function locates and returns
        the GeneID if present, or returns `None` otherwise.
        """
        values = self.get_attribute('Dbxref', as_list=True)
        if values is None:
            return None
        for value in values:
            if value.startswith('GeneID:'):
                key, geneid = value.split(':')
                return geneid
        return None

    @property
    def cdslen(self):
        """
        Translated length of this feature.

        Undefined for non-mRNA features.
        """
        if self.type != 'mRNA':
            return None

        return sum([len(c) for c in self.children if c.type == 'CDS'])

    def overlap(self, rng):
        """Report whether this feature overlaps with the specified range."""
        return self._range.overlap(rng)

    def contains(self, rng):
        """Report whether this feature contains the specified range."""
        return self._range.contains(rng)

    def contains_point(self, point):
        """Report whether this feature contains the specified point."""
        return self._range.contains_point(point)
