
COM = fvc::domainIntegrate(mesh.C() * alpha1)/fvc::domainIntegrate(alpha1);
KE = 0.5 * fvc::domainIntegrate(rho * (U & U));

SpuriousFile << runTime.timeName() << '\t' << COM.value().x() << '\t' << COM.value().y() << '\t' << COM.value().z() << '\t' << KE.value() << nl;
