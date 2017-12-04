
const scalar alphaTol(dict.lookupOrDefault<scalar>("alphaTol", 1e-8));
phaseidentifier = advector.phasedetermination();

forAll(alpha1,celli){
    if(alpha1[celli] > 1 - alphaTol)
	liquidorvapor[celli] = -1;
    else
    if(alpha1[celli] < alphaTol)
	liquidorvapor[celli] = 1;
    else
	if(phaseidentifier[celli])
	    liquidorvapor[celli] = -1;
	else
	    liquidorvapor[celli] = 1;
}
