function output= generateJVectors(path, sliceT, pm, pn, sst, Qo, EP, SW, rho0, Cp, tx, ty,...
                                tmag, dxf, dyf, theta_s, theta_b, hc, h,zl, g, f, zmin, OMEGAZ)
alpha = 2e-4;
beta  = -7e-5;
thmix = GetVarROMS(path, 0, {'temp_hmix', '(1)'}, sliceT);
tvmix = GetVarROMS(path, 0, {'temp_vmix', '(1)'}, sliceT);
shmix = GetVarROMS(path, 0, {'salt_hmix', '(1)'}, sliceT);
svmix = GetVarROMS(path, 0, {'salt_vmix', '(1)'}, sliceT);
D     = g*alpha.*rho0.*(thmix + tvmix)   + g*beta.*rho0.*(shmix + svmix);
output.Jd = OMEGAZ.*D;

end