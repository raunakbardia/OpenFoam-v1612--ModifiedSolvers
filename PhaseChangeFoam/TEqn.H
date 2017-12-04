
rhoCv = alpha1 * rho1 * Cvl + (1 - alpha1) * rho2 * Cvv;

rhoCp = alpha1 * rho1 * Cpl + (1 - alpha1) * rho2 * Cpv;

Hphi = fvc::interpolate(rhoCp) * phi;

fvScalarMatrix TEqn
(
    fvm::ddt(rhoCv,T)
    + fvm::div(Hphi,T)
    - fvm::laplacian(k,T)
);
TEqn.solve();
~                                                                                                                                      
~                                                                                                                                      
~                                                                                                                                      
~                                                                                                                                      
~   
