NEURON {
	SUFFIX cdiffdemo
	USEION A WRITE Ai VALENCE 2
	RANGE molA, rand
}

UNITS {
    (mM) = (milli/liter)
    (um) = (micron)
    FARADAY = (faraday) (coulomb)
    PI = (pi) (1)
}


ASSIGNED {
	: current concentration of A in encapsulating tet
	molA
    : random number to determine section output
    rand
}


STATE {
	Ai (mM)
}

INITIAL {
    molA = 0.0
    Ai = molA
}

BREAKPOINT {
    Ai = molA*rand
}