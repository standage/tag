#!/usr/bin/env python
#
# -----------------------------------------------------------------------------
# Copyright (C) 2015 Daniel Standage <daniel.standage@gmail.com>
#
# This file is part of aeneas (http://github.com/standage/aeneas) and is
# licensed under the ISC license: see LICENSE.txt.
# -----------------------------------------------------------------------------
"""Represents a feature entry from a GFF3 file."""

from .region import Region


class Feature(object):

    def __init__(self, data):
        fields = data.split('\t')
        assert len(fields) == 9
        self.children = None
        self.multi_rep = None
        self.siblings = None

        self._seqid = fields[0]
        self._source = fields[1]
        self._type = fields[2]
        self._region = Region(int(fields[3]), int(fields[4]))
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

        assert self._strand in ['+', '-', '.'], 'invalid strand "%s"' % \
            self._strand
        if self.phase is not None:
            assert self.phase in [0, 1, 2], 'invalid phase "%s"' % self.phase

    def __str__(self):
        """String representation of the feature, sans children."""
        score = '.'
        if self.score is not None:
            score = "%.3f" % self.score
        phase = '.'
        if self.phase is not None:
            phase = '%d' % self.phase
        return '\t'.join([
            self.seqid, self.source, self.type, '%d' % self.start,
            '%d' % self.end, score, self.strand, phase, self.attributes
        ])

    def __repr__(self):
        """Full representation of the feature, with children"""
        string = ''
        for feature in self:
            if string != '':
                string += '\n'
            string += feature.__str__()
        return string

    def __len__(self):
        return len(self._region)

    def __lt__(self, other):
        """Rich comparison operator for Python 3 support."""
        if self.seqid < other.seqid:
            return True
        elif self.seqid > other.seqid:
            return False
        return self._region.__lt__(other._region)

    def __le__(self, other):
        """Rich comparison operator for Python 3 support."""
        if self.seqid < other.seqid:
            return True
        elif self.seqid > other.seqid:
            return False
        else:
            return self._region.__le__(other._region)

    def __gt__(self, other):
        """Rich comparison operator for Python 3 support."""
        if self.seqid > other.seqid:
            return True
        elif self.seqid < other.seqid:
            return False
        else:
            return self._region.__gt__(other._region)

    def __ge__(self, other):
        """Rich comparison operator for Python 3 support."""
        if self.seqid > other.seqid:
            return True
        elif self.seqid < other.seqid:
            return False
        else:
            return self._region.__ge__(other._region)

    def __iter__(self):
        """
        Perform a depth-first search of the feature graph.

        Nodes with multiple parents are reported only once.
        """
        yield self
        if self.children is not None:
            already_seen = dict()
            for child in self.children:
                for subchild in child:
                    if subchild not in already_seen:
                        already_seen[subchild] = 1
                        yield subchild

    def add_child(self, child, regioncheck=True):
        if regioncheck is True:
            assert self.seqid == child.seqid and self._region.contains(child)
            assert self._strand == child._strand
        if self.children is None:
            self.children = list()
        self.children.append(child)
        self.children.sort()

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
        assert start > 0 and end > 0 and start <= end
        self._region = Region(start, end)

    def transform(self, offset, newseqid=None):
        for feature in self:
            feature._region.transform(offset)
            if newseqid is not None:
                feature._seqid = newseqid

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
            attrs.append('%s=%s' % (attrkey, value))
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
                    if oldid is not None:
                        child.add_attribute('Parent', attrvalue,
                                            oldvalue=oldid)
                    else:
                        child.add_attribute('Parent', attrvalue)
            self._attrs[attrkey] = attrvalue
            return

        # Handle all other attribute types
        if oldvalue is not None:
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
            key, value = kvp.split('=')
            if key == 'ID':
                assert ',' not in value
                attributes[key] = value
                continue
            values = value.split(',')
            valdict = dict((val, True) for val in values)
            attributes[key] = valdict
        return attributes


# -----------------------------------------------------------------------------
# Unit tests
# -----------------------------------------------------------------------------

def eden():
    """
    Data fixture for unit tests.

    This function is not a good example of how to parse arbitrary GFF3 files,
    it's quick and dirty for this particular case.
    """
    edenfile = open('testdata/eden.gff3', 'r')
    edengff3 = edenfile.readlines()
    gene = None
    mRNAs = list()
    cdsreps = dict()
    for line in edengff3:
        if line.startswith('#'):
            continue
        feature = Feature(line.rstrip())
        if feature.type == 'gene':
            gene = feature
        elif feature.type == 'TF_binding_site':
            gene.add_child(feature)
        elif feature.type == 'mRNA':
            mRNAs.append(feature)
            gene.add_child(feature)
        else:
            for parent in feature.get_attribute('Parent', as_list=True):
                for mRNA in mRNAs:
                    if mRNA.get_attribute('ID') == parent:
                        mRNA.add_child(feature)
                        break
            if feature.type == 'CDS':
                cdsid = feature.get_attribute('ID')
                if cdsid not in cdsreps:
                    cdsreps[cdsid] = feature
                else:
                    cdsreps[cdsid].add_sibling(feature)
    assert len(gene.children) == 4
    return gene


def test_basic():
    """[aeneas::Feature] Test basic constructor."""
    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', '.']
    f1 = Feature('\t'.join(gff3))
    assert '%s' % f1 == '\t'.join(gff3)
    assert len(f1) == 1001

    gff3[8] = 'ID=gene1;Name=EDEN'
    f2 = Feature('\t'.join(gff3))
    assert '%s' % f2 == '\t'.join(gff3)

    gff3[5] = '0.28466'
    f3 = Feature('\t'.join(gff3))
    assert abs(f3.score - 0.28466) <= 0.00001
    gff3[5] = '0.285'
    assert '%s' % f3 == '\t'.join(gff3)


def test_seqid():
    """[aeneas::Feature] Test seqid handling."""
    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', 'ID=g1']
    f1 = Feature('\t'.join(gff3))
    assert f1.seqid == 'chr'

    f1.seqid = 'scaffold_28'
    assert f1.seqid == 'scaffold_28'
    gff3[0] = 'scaffold_28'
    assert '%s' % f1 == '\t'.join(gff3)

    gff3 = ['scaffold_28', 'vim', 'mRNA', '1000', '2000', '.', '+', '.',
            'ID=t1;Parent=g1']
    f2 = Feature('\t'.join(gff3))
    gff3 = ['scaffold_28', 'vim', 'exon', '1000', '2000', '.', '+', '.',
            'Parent=t1']
    f3 = Feature('\t'.join(gff3))
    f2.add_child(f3)
    f1.add_child(f2)
    for feature in f1:
        assert feature.seqid == 'scaffold_28'
    f1.seqid = 'ctg123'
    for feature in f1:
        assert feature.seqid == 'ctg123'


def test_source():
    """[aeneas::Feature] Test source handling."""
    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', 'ID=g1']
    f1 = Feature('\t'.join(gff3))
    assert f1.source == 'vim'

    f1.source = 'nano'
    assert f1.source == 'nano'
    gff3[1] = 'nano'
    assert '%s' % f1 == '\t'.join(gff3)

    gff3 = ['chr', 'emacs', 'mRNA', '1000', '2000', '.', '+', '.',
            'ID=t1;Parent=g1']
    f2 = Feature('\t'.join(gff3))
    gff3 = ['chr', 'nano', 'exon', '1000', '2000', '.', '+', '.',
            'Parent=t1']
    f3 = Feature('\t'.join(gff3))
    f2.add_child(f3)
    f1.add_child(f2)
    for feature in f1:
        if feature.type == 'mRNA':
            assert feature.source == 'emacs'
        else:
            assert feature.source == 'nano'
    f1.source = 'gedit'
    for feature in f1:
        if feature.type == 'mRNA':
            assert feature.source == 'emacs'
        else:
            assert feature.source == 'gedit', '%s=%s' % (feature.type,
                                                         feature.source)


def test_type():
    """[aeneas::Feature] Test type handling."""
    gff3 = ['chr', 'vim', 'mRNA', '1000', '2000', '.', '+', '.', 'ID=mRNA1']
    f1 = Feature('\t'.join(gff3))
    assert f1.type == 'mRNA'

    gff3[2] = 'messenger RNA'
    f2 = Feature('\t'.join(gff3))
    assert f2.type == 'messenger RNA'


def test_region():
    """[aeneas::Feature] Test coordinate handling."""
    gff3 = ['chr', 'vim', 'gene', '1000', '2000', '.', '+', '.', 'ID=g1']
    f1 = Feature('\t'.join(gff3))
    assert f1._region == Region(1000, 2000)
    assert f1.start == 1000 and f1.end == 2000

    gff3 = ['contig5', 'vim', 'mRNA', '500', '2500', '.', '+', '.',
            'ID=t1;Parent=g1']
    f2 = Feature('\t'.join(gff3))
    try:
        f1.add_child(f2)
    except AssertionError:
        pass
    f1.add_child(f2, regioncheck=False)
    assert len(f1.children) == 1
    assert f2._region == Region(500, 2500)

    f2.set_coord(1000, 2000)
    assert f2.start == 1000 and f2.end == 2000

    f1.transform(100000)
    for feature in f1:
        assert feature._region == Region(101000, 102000), \
            '%s %r' % (feature.type, feature._region)

    f1.transform(100000, newseqid='scf89')
    for feature in f1:
        assert feature.seqid == 'scf89'
        assert feature._region == Region(201000, 202000), \
            '%s %r' % (feature.type, feature._region)


def test_score():
    """[aeneas::Feature] Test score handling."""
    gff3 = ['chr', 'vim', 'EST_match', '57229', '57404', '.', '+', '.', '.']
    f1 = Feature('\t'.join(gff3))
    assert f1.score is None

    gff3[5] = '0.97'
    f2 = Feature('\t'.join(gff3))
    assert abs(f2.score - 0.97) <= 0.00001
    gff3[5] = '0.970'
    assert '%s' % f2 == '\t'.join(gff3)

    gff3[5] = '-1.8332'
    f3 = Feature('\t'.join(gff3))
    assert abs(f3.score + 1.8332) <= 0.00001
    gff3[5] = '-1.833'
    assert '%s' % f3 == '\t'.join(gff3)


def test_strand():
    """[aeneas::Feature] Test strand handling."""
    gff3 = ['chr', 'vim', 'EST_match', '57229', '57404', '.', '.', '.', '.']
    f1 = Feature('\t'.join(gff3))
    assert f1.strand == '.'
    assert '%s' % f1 == '\t'.join(gff3)

    gff3[6] = '-'
    f2 = Feature('\t'.join(gff3))
    assert f2.strand == '-'
    assert '%s' % f2 == '\t'.join(gff3)

    gff3[2:5] = 'match_part', '57229', '57298'
    f3 = Feature('\t'.join(gff3))
    assert f2.strand == '-'

    gff3[3:5] = '57377', '57404'
    gff3[6] = '+'
    f4 = Feature('\t'.join(gff3))
    assert f4.strand == '+'

    f2.add_child(f3)
    try:
        f2.add_child(f4)
    except AssertionError:
        pass
    assert len(f2.children) == 1
    f4._strand = '-'
    f2.add_child(f4)
    assert len(f2.children) == 2

    gff3 = ['chr', 'vim', 'EST_match', '57229', '57404', '.', '~', '.', '.']
    try:
        f5 = Feature('\t'.join(gff3))
    except AssertionError:
        pass


def test_phase():
    """[aeneas::Feature] Test phase handling."""
    gff3 = ['chr', 'vim', 'mRNA', '1001', '1420', '.', '+', '.', 'ID=t1']
    m1 = Feature('\t'.join(gff3))
    assert m1.phase is None

    gff3 = ['chr', 'vim', 'CDS', '1001', '1100', '.', '+', '0', 'ID=CDS1']
    c1 = Feature('\t'.join(gff3))
    assert c1.phase == 0

    gff3[3:5] = '1201', '1236'
    gff3[7] = '2'
    c2 = Feature('\t'.join(gff3))
    c1.add_sibling(c2)
    assert c2.phase == 2

    gff3[3:5] = '1301', '1364'
    c3 = Feature('\t'.join(gff3))
    c1.add_sibling(c3)
    assert c3.phase == 2

    gff3[3:5] = '1401', '1420'
    gff3[7] = '1'
    c4 = Feature('\t'.join(gff3))
    c1.add_sibling(c4)
    assert c4.phase == 1


def test_attributes():
    """[aeneas::Feature] Test attribute handling."""
    gene = eden()
    assert gene.get_attribute('ID') == 'gene00001'
    assert gene.get_attribute('ID', as_list=True) == 'gene00001'
    assert gene.get_attribute('Name') == 'EDEN'
    assert gene.get_attribute('Name', as_list=True) == ['EDEN']
    assert gene.get_attribute_keys() == ['ID', 'Name']

    gene.add_attribute('Name', 'Gandalf')
    assert gene.get_attribute('Name') == 'Gandalf'
    gene.add_attribute('Name', 'Aragorn', append=True)
    assert gene.get_attribute('Name') == ['Aragorn', 'Gandalf']
    assert gene.get_attribute('Name', as_string=True) == 'Aragorn,Gandalf'

    gene.add_attribute('ID', 'g1')
    for child in gene.children:
        assert child.get_attribute('Parent') == 'g1'

    for feature in gene:
        if feature.type == 'exon' and feature._region == Region(5000, 5500):
            assert feature.get_attribute('Parent') == [
                'mRNA00001',
                'mRNA00002',
                'mRNA00003']
            assert feature.get_attribute('Parent', as_string=True) == \
                'mRNA00001,mRNA00002,mRNA00003'

    for feature in gene:
        if feature.get_attribute('ID') == 'mRNA00003':
            feature.add_attribute('ID', 'mRNA3')
        elif feature.get_attribute('ID') in ['cds00003', 'cds00004']:
            assert feature.get_attribute('Parent') == 'mRNA3'

    assert '%r' % gene == open('testdata/eden-mod.gff3', 'r').read().rstrip()

    gene.add_attribute('Note', 'I need to test access of other attributes')
    assert gene.attributes == ('ID=g1;Name=Aragorn,Gandalf;Note='
                               'I need to test access of other attributes')

    gff3 = ['chr', 'vim', 'mRNA', '1001', '1420', '.', '+', '.', '.']
    m1 = Feature('\t'.join(gff3))
    assert m1.phase is None

    gff3 = ['chr', 'vim', 'CDS', '1001', '1420', '.', '+', '.', '.']
    c1 = Feature('\t'.join(gff3))
    m1.add_child(c1)
    m1.add_attribute('ID', 'mRNA1')
    assert c1.get_attribute('Parent') == 'mRNA1'


def test_multi():
    """
    [aeneas::Feature] Test handling of multi-features.

    In GFF3 discontiguous sequence features are often encoded using multi-
    features: that is, a single feature described across mulitple entries
    (lines). This unit test validates the handling of multi-features.
    """
    gene = eden()
    for feature in gene:
        if feature.type == 'CDS':
            assert feature.is_multi
            if feature.get_attribute('ID') == 'cds00001':
                assert feature.multi_rep._region == Region(1201, 1500)
            elif feature.get_attribute('ID') == 'cds00002':
                assert feature.multi_rep._region == Region(1201, 1500)
            elif feature.get_attribute('ID') == 'cds00003':
                assert feature.multi_rep._region == Region(3301, 3902)
            elif feature.get_attribute('ID') == 'cds00004':
                assert feature.multi_rep._region == Region(3391, 3902)
        else:
            assert not feature.is_multi and feature.multi_rep is None

    for feature in gene:
        if feature.get_attribute('ID') == 'cds00004':
            if feature.multi_rep == feature:
                feature.type = 'coding sequence'
            assert feature.type == 'coding sequence', feature._region


def test_compare():
    """[aeneas::Feature] Test comparison and sorting of features."""
    gff3 = ['chr1', 'vim', 'gene', '1000', '2000', '.', '+', '.', '.']
    g1 = Feature('\t'.join(gff3))
    gff3 = ['chr1', 'vim', 'gene', '3000', '4000', '.', '+', '.', '.']
    g2 = Feature('\t'.join(gff3))
    gff3 = ['chr1', 'vim', 'gene', '3000', '4500', '.', '+', '.', '.']
    g3 = Feature('\t'.join(gff3))

    assert g1.__lt__(g2) is True
    assert g1.__le__(g3) is True
    assert g3.__gt__(g1) is True
    assert g1.__gt__(g3) is False
    assert g3.__ge__(g2) is True
    assert g3.__le__(g2) is False
    assert sorted([g3, g2, g1]) == [g1, g2, g3]

    gff3 = ['chr10', 'vim', 'gene', '100', '400', '.', '-', '.', '.']
    g4 = Feature('\t'.join(gff3))
    gff3 = ['chr2', 'vim', 'gene', '2000', '2500', '.', '-', '.', '.']
    g5 = Feature('\t'.join(gff3))

    assert g2.__le__(g4) is True
    assert g4.__le__(g2) is False
    assert g4.__lt__(g5) is True
    assert g5.__lt__(g4) is False
    assert g5.__gt__(g1) is True
    assert g1.__gt__(g5) is False
    assert g5.__ge__(g1) is True
    assert g1.__ge__(g5) is False
    assert g5.__ge__(g5) is True
    assert g5.__le__(g5) is True
    assert sorted([g3, g5, g1, g4, g2]) == [g1, g2, g3, g4, g5]
