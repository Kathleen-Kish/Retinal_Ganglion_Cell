TITLE decay of submembrane calcium concentration
:
: Internal calcium concentration due to calcium currents and pump.
: Pump action is approximated by a simple first order decay.

NEURON {
	SUFFIX cad
	USEION ca READ ica, cai WRITE cai
	RANGE depth,kd,cainf,taur
}

UNITS {
	(molar) = (1/liter)			: moles do not appear in units
	(mM)	= (millimolar)
	(um)	= (micron)
	(mA)	= (milliamp)
	(msM)	= (ms mM)
}

CONSTANT {
	FARADAY = 96489		(coul)		: moles do not appear in units
}

PARAMETER {
	depth	= .1	(um)		: depth of shell
	taur	= 1.5	(ms)		: remove first-order decay
	cainf	= 0.0001 (mM)
	kd	= 0.0001	(mM)
}

STATE {
	cai		(mM) 
}

INITIAL {
	cai = kd
}

ASSIGNED {
	ica		(mA/cm2)
	drive_channel	(mM/ms)
}
	
BREAKPOINT {
	SOLVE state METHOD derivimplicit
}

DERIVATIVE state { 

	drive_channel =  - (10000) * ica / (2 * FARADAY * depth)

	if (drive_channel <= 0.) { drive_channel = 0. }	: cannot pump below resting level

	cai' = drive_channel + (cainf-cai)/taur
}

