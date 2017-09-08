COM_y = fvc::domainIntegrate(mesh.C() * alpha1)/fvc::domainIntegrate(alpha1);

COMFile	<< runTime.timeName() << '\t' << COM_y.value().y() << nl;
