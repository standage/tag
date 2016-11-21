#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the BSD 3-clause license: see LICENSE.txt.
# -----------------------------------------------------------------------------

from .comment import Comment
from .directive import Directive
from .range import Range
from .sequence import Sequence


class Feature(object):
    """Represents a feature entry from a GFF3 file."""

    def __init__(self, data):
        fields = data.split('\t')
        assert len(fields) == 9
        self.children = None
        self.multi_rep = None
        self.siblings = None

        self._seqid = fields[0]
        self._source = fields[1]
        self._type = fields[2]
        self._region = Range(int(fields[3]) - 1, int(fields[4]))
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
        if self.children is not None:
            string += '\n###'
        return string

    def __len__(self):
        return len(self._region)

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
        elif self._region.__eq__(other._region):
            return self.type > other.type
        return self._region.__lt__(other._region)

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
        elif self._region.__eq__(other._region):
            return self.type >= other.type
        return self._region.__le__(other._region)

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
        elif self._region.__eq__(other._region):
            return self.type < other.type
        return self._region.__gt__(other._region)

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
        elif self._region.__eq__(other._region):
            return self.type <= other.type
        return self._region.__ge__(other._region)

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
            if self.children is not None:
                for child in reversed(self.children):
                    child._visit(L, marked, tempmarked)
            marked[self] = True
            del tempmarked[self]
            L.insert(0, self)

    def add_child(self, child, rangecheck=False):
        assert self.seqid == child.seqid, \
            (
                'seqid mismatch for feature {} ({} vs {})'.format(
                    self.fid, self.seqid, child.seqid
                )
            )
        if rangecheck is True:
            assert self._strand == child._strand, \
                ('child of feature {} has a different strand'.format(self.fid))
            assert self._region.contains(child._region), \
                (
                    'child of feature {} is not contained within its span '
                    '({}-{})'.format(self.fid, child.start, child.end)
                )
        if self.children is None:
            self.children = list()
        self.children.append(child)
        self.children.sort()

    @property
    def fid(self):
        return self.get_attribute('ID')

    @property
    def is_multi(self):
        return self.multi_rep is not None

    def add_sibling(self, sibling):
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
        return self._region.start

    @property
    def end(self):
        return self._region.end

    def set_coord(self, start, end):
        self._region = Range(start, end)

    def transform(self, offset, newseqid=None):
        for feature in self:
            feature._region.transform(offset)
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
        Attributes stored as nested dictionaries.

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

    def get_attribute_keys(self):
        return sorted(list(self._attrs))

    def parse_attributes(self, attrstring):
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
