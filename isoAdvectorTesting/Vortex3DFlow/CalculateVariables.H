const scalar period(3);
const scalar t(runTime.time().value());
const scalar dt(runTime.deltaT().value());
const scalar uFactor
(
	.5*(Foam::cos(M_PI*t/period)
        + Foam::cos(M_PI*(t + dt)/period))
);
phi = uFactor*phi0;
U = uFactor*U0;
