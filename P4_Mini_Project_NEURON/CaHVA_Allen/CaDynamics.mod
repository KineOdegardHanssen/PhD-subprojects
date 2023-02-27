: Dynamics that track inside calcium concentration
: modified from Destexhe et al. 1994

NEURON	{
	SUFFIX CaDynamics
	USEION ca READ ica WRITE cai, cao
	RANGE decay, gamma, minCai, depth
}

UNITS	{
	(mV) = (millivolt)
	(mA) = (milliamp)
	FARADAY = (faraday) (coulombs)
	(molar) = (1/liter)
	(mM) = (millimolar)
	(um)	= (micron)
}

PARAMETER	{
	gamma = 1.0 : percent of free calcium (not buffered) 0.05 :
	decay = 80 (ms) : rate of removal of calcium
	depth = 0.1 (um) : depth of shell
	minCai = 1e-4 (mM)
        cao0 = 2 (mM)
}

ASSIGNED	{ica (mA/cm2)}

INITIAL {
	cai = minCai
        cao = cao0
}

STATE	{
	cai (mM)
        cao (mM)
}

BREAKPOINT	{ SOLVE states METHOD cnexp }

DERIVATIVE states	{
	cai' = -(10000)*(ica*0.2/(2*FARADAY*depth)) - (cai - minCai)/5.0
        cao' = 0
}
