########################################################################
# USED TO SPECIFY COLUMN LABELS...REPLACE WITH REGEXP PATTERNS,
# MOVE TO SEQUENCES MODULE.
########################################################################
SITES = [
	"Site_P4", "Site_P3", "Site_P2", "Site_P1", "Site_P1prime",
	"Site_P2prime", "Site_P3prime", "Site_P4prime",
]

MEROPS_POSITIONS = [
	"p4","p3","p2","p1","p1'","p2'","p3'","p4'",
]

PHOSPHO_POSITIONS = [
	"p1","p2","p3","p4","p5","p6","p7",
	"p8","p9","p10","p11","p12","p13",
]

POSITIONS = MEROPS_POSITIONS
IGNORE_POSITIONS = ['p7']
########################################################################

THREE_LETTER_CODES = [
    'Ala','Arg','Asn','Asp','Asx','Cys','Glu','Gln','Glx',
    'Gly','His','Ile','Leu','Lys','Met','Phe','Pro','Ser',
    'Thr','Trp','Tyr','Val',
]

########################################################################
# MAPS SINGLE-LETTER CATALYTIC TYPE CODES FROM MEROPS IDs TO DESCRIPTIVE
# CATALYTIC TYPE NAMES. ONLY USED FOR GENERATION OF LABELS ON PLOTS.
########################################################################
CATALYTIC_TYPES = {
	'A':'Aspartic', 'C':'Cysteine','G':'Glutamic','M':'Metallo',
	'N':'Asparagine','P':'Mixed','S':'Serine','T':'Threonine',
	'U':'Unknown',
}
########################################################################