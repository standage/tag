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


class Feature(object):
    """
    Represents a feature entry from a GFF3 file.

    >>> feature = tag.demo_feature()
    >>> feature.seqid
    'contig1'
    >>> feature.source
    'snap'
    >>> feature.type
    'gene'
    >>> feature.start, feature.end
    (999, 7500)
    >>> feature.score is None
    True
    >>> feature.strand
    '+'
    >>> feature.phase is None
    True
    >>> feature.attributes
    'ID=gene1'
    >>> feature.num_children
    1
    >>> feature.is_multi
    False
    >>> feature.is_toplevel
    True
    >>> for child in feature:
    ...     if child.type == 'CDS':
    ...         assert child.get_attribute('ID') == 'cds1'
    >>> feature.slug
    'gene@contig1[1000, 7500]'
    """

    def __init__(self, data):
        fields = data.split('\t')
        assert len(fields) == 9
        self.children = None
        self.multi_rep = None
        self.siblings = None

        self._seqid = fields[0]
        self._source = fields[1]
        self._type = fields[2]
        self._range = Range(int(fields[3]) - 1, int(fields[4]))
        if fields[5] == '.':
            self.score = None
        else:
            self.score = float(fields[5])
        self._strand = fields[6]
        if fields[7] == '.':
            self.phase = None
        else:
            self.phase = int(fields[7])
        self._attrs = self.parse_attributes(fields[8])

        assert self._strand in ['+', '-', '.'], \
            'invalid strand "{}"'.format(self._strand)
        if self.phase is not None:
            assert self.phase in [0, 1, 2], \
                'invalid phase "{}"'.format(self.phase)

    def __str__(self):
        """String representation of the feature, sans children."""
        score = '.'
        if self.score is not None:
            score = "{:.3f}".format(self.score)
        phase = '.'
        if self.phase is not None:
            phase = str(self.phase)
        return '\t'.join([
            self.seqid, self.source, self.type, str(self.start + 1),
            str(self.end), score, self.strand, phase, self.attributes
        ])

    def __repr__(self):
        """Full representation of the feature, with children"""
        string = ''
        for feature in self:
            if string != '':
                string += '\n'
            string += str(feature)
        if self.children is not None or self.is_multi:
            string += '\n###'
        return string

    def __len__(self):
        return len(self._range)

    def __lt__(self, other):
        if isinstance(other, Directive) or isinstance(other, Comment):
            return False
        elif isinstance(other, Sequence):
            return True
        assert isinstance(other, Feature)

        if self.seqid < other.seqid:
            return True
        elif self.seqid > other.seqid:
            return False
        elif self._range == other._range:
            return self.type > other.type
        return self._range < other._range

    def __le__(self, other):
        if isinstance(other, Directive) or isinstance(other, Comment):
            return False
        elif isinstance(other, Sequence):
            return True
        assert isinstance(other, Feature)

        if self.seqid < other.seqid:
            return True
        elif self.seqid > other.seqid:
            return False
        elif self._range == other._range:
            return self.type >= other.type
        return self._range <= other._range

    def __gt__(self, other):
        if isinstance(other, Directive) or isinstance(other, Comment):
            return True
        elif isinstance(other, Sequence):
            return False
        assert isinstance(other, Feature)

        if self.seqid > other.seqid:
            return True
        elif self.seqid < other.seqid:
            return False
        elif self._range == other._range:
            return self.type < other.type
        return self._range > other._range

    def __ge__(self, other):
        if isinstance(other, Directive) or isinstance(other, Comment):
            return True
        elif isinstance(other, Sequence):
            return False
        assert isinstance(other, Feature)

        if self.seqid > other.seqid:
            return True
        elif self.seqid < other.seqid:
            return False
        elif self._range == other._range:
            return self.type <= other.type
        return self._range >= other._range

    def __iter__(self):
        """Generator iterates through a feature and all its subfeatures."""
        sorted_features = list()
        self._visit(L=sorted_features, marked={}, tempmarked={})
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
        if self in tempmarked:
            raise Exception('feature graph is cyclic')
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
    def is_toplevel(self):
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
    def strand(self):
        return self._strand

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
        if attrstring == '.':
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
