intL = fvc::domainIntegrate(mag(alpha1 - alpha0));
IntL << runTime.timeName() << '\t' << intL.value() << nl;
