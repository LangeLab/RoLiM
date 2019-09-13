from collections import namedtuple

CompoundResidue = namedtuple('CompoundResidue',
								['code', 'residues', 'vector', 'frequency'])

compound_residues = {
	'helix breaker': CompoundResidue(
		code='1',
		residues=[
			'G',
			'P',
		],
		vector=None,
		frequency=None
	),
	'tiny (flexible)': CompoundResidue(
		code='2',
		residues=[
			'G',
		]
	),
	'no nh': CompoundResidue(
		code='3',
		residues=[
			'P',
		]
	),
	'aliphatic': CompoundResidue(
		code='4',
		residues=[
			'P',
			'A',
			'L',
			'V',
			'I',
		]
	),
	'beta-branched (rigid)': CompoundResidue(
		code='5',
		residues=[
			'V',
			'I',
			'T',
		]
	),
	'hydrophobic': CompoundResidue(
		code='6',
		residues=[
			'M',
			'F',
			'W',
			'Y',
			'P', 
			'A',
			'L',
			'V',
			'I',
			'T',
			'C',
		]
	),
	'disulfide': CompoundResidue(
		code='7',
		residues=[
			'C',
		]
	),
	'nucleophile in proteases': CompoundResidue(
		code='8',
		residues=[
			'C',
			'S',
		]
	),
	'polar': CompoundResidue(
		code='9',
		residues=[
			'N',
			'Q',
			'S',
			'C',
			'T',
			'Y',
			'H',
			'R',
			'K',
			'D',
			'E',
		]
	),
	'aromatic': CompoundResidue(
		code='10',
		residues=[
			'F',
			'W',
			'Y',
			'H',
		]
		
	),
	'proton donor/acceptor in catalysis, Fe2+, Zn2+ binder': CompoundResidue(
		code='11',
		residues=[
			'C',
			'H',
		]
	),
	'base, DNA, RNA binder': CompoundResidue(
		code='12',
		residues=[
			'R',
			'K',
		]
	),
	'acid, Ca2+, Mg2+ binder': CompoundResidue(
		code='13',
		residues=[
			'D',
			'E',
		]

	),
	'charged': CompoundResidue(
		code='14',
		residues=[
			'H',
			'R',
			'K',
			'D',
			'E',
		]
	),
}