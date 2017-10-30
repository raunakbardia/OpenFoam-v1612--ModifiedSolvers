KE = 0.5 * fvc::domainIntegrate(rho * (U & U));

KEFile << runTime.timeName() << '\t' << KE.value() << nl;
