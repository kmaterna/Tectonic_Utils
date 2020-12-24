
# Calculations that convert between moment and magnitude, etc. 

import numpy as np 


def moment_from_muad(mu, A, d):
	# moment = mu * A * d
	return mu*A*d;

def mw_from_moment(moment):
	# Takes newton meters, returns a moment magnitude
	moment = moment*1e7;
	mw = (2/3)*np.log10(moment) - 10.7
	return mw;

def moment_from_mw(Mw):
	# Current definition of moment magnitude, returning moment in newton meters
	exponent = 1.5*Mw + 1.5*10.7;
	moment = np.power(10, exponent);
	moment_newton_meters = moment * 1e-7;
	return moment_newton_meters;
