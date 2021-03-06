function D = bRHS(TRHS, SRHS, T, S, rho0, g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% BASED ON rho_eos from CROCO ROMS tools
% %  Further Information:  
%  http://www.brest.ird.fr/Roms_tools/
%
% 
% Modified by J. Wenegrat
% Calculates a derivative operator applied to the nonlinear buoyancy
% equation in terms of the same operator applied to T and S.
%
% ex. db/dt  can be calculated as bRHS(dTdt, dSdt, T, S, rho0, g)
% Simple to show that the RHS of the buoyancy equation is equivalent to
% bRHS(TRHS, SRHS, T, S, rho0, g)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% z_r=0;
% A00=+19092.56;A01=+209.8925;
% A02=-3.041638;A03=-1.852732e-3;A04=-1.361629e-5;A10=104.4077;
% A11=-6.500517;A12=+0.1553190;A13=2.326469e-4;AS0=-5.587545;
% AS1=+0.7390729;AS2=-1.909078e-2;B00=+4.721788e-1;B01=+1.028859e-2;
% B02=-2.512549e-4;B03=-5.939910e-7;B10=-1.571896e-2;B11=-2.598241e-4;
% B12=+7.267926e-6;BS1=+2.042967e-3;E00=+1.045941e-5;E01=-5.782165e-10;
% E02=+1.296821e-7;E10=-2.595994e-7;E11=-1.248266e-9;E12=-3.508914e-9;

QR=+999.842594;Q01=+6.793952e-2;Q02=-9.095290e-3;
Q03=+1.001685e-4;Q04=-1.120083e-6;Q05=+6.536332e-9;Q10=+0.824493;
Q11=-4.08990e-3;Q12=+7.64380e-5;Q13=-8.24670e-7;Q14=+5.38750e-9;
QS0=-5.72466e-3;QS1=+1.02270e-4;QS2=-1.65460e-6;Q20=+4.8314e-4;


% sqrtTs=sqrt(Ts);
% 
% K0=A00+Tt.*(A01+Tt.*(A02+Tt.*(A03+Tt.*A04)))...
%    +Ts.*(A10+Tt.*(A11+Tt.*(A12+Tt.*A13))...
%    +sqrtTs.*(AS0+Tt.*(AS1+Tt.*AS2)));
% K1=B00+Tt.*(B01+Tt.*(B02+Tt.*B03))...
%    +Ts.*(B10+Tt.*(B11+Tt.*B12)+sqrtTs.*BS1);
% K2=E00+Tt.*(E01+Tt.*E02)...
%    +Ts.*(E10+Tt.*(E11+Tt.*E12));
% rho1=QR+Tt.*(Q01+Tt.*(Q02+Tt.*(Q03+Tt.*(Q04+Tt.*Q05))))...
%      +Ts.*(Q10+Tt.*(Q11+Tt.*(Q12+Tt.*(Q13+Tt.*Q14)))...
%      +sqrtTs.*(QS0+Tt.*(QS1+Tt.*QS2))+Ts.*Q20);
% rho=rho1./(1+0.1.*z_r./(K0-z_r.*(K1-z_r.*K2)));


% 
% T1 = Tt.*Q01
% T2 = Tt.^2.*Q02
% T3 = Tt.^3.*Q03
% T4 = Tt.^4.*Q04
% T5 = Tt.^5.*Q05
% 
% S1 = Ts.*Q10 
% S2 = sqrt(Ts).*Ts.*QS0
% S3 = sqrt(Ts).*Ts.*Tt.*QS1
% S4 = sqrt(Ts).*Ts.*Tt.^2.*QS2
% S5 = Ts.^2*Q20
% 
% M1 = Ts.*Tt.*Q11
% M2 = Ts.*Tt.^2.*Q12
% M3 = Ts.*Tt.^3.*Q13
% M4 = Ts.*Tt.^4.*Q14



T1p = TRHS.*Q01;
T2p = 2.*T.*TRHS.*Q02;
T3p = 3.*T.^2.*TRHS.*Q03;
T4p = 4.*T.^3.*TRHS.*Q04;
T5p = 5.*T.^4.*TRHS.*Q05;

S1p = SRHS.*Q10;
S2p = 3/2.*sqrt(S).*SRHS.*QS0;
S3p = (3/2*sqrt(S).*T.*SRHS + sqrt(S).*S.*TRHS).*QS1;
S4p = (3/2.*sqrt(S).*T.^2.*SRHS + 2*sqrt(S).*T.*S.*TRHS).*QS2;
S5p = 2*S.*SRHS.*Q20;

M1p = (T.*SRHS + S.*TRHS).*Q11;
M2p = (T.^2.*SRHS + 2*T.*S.*TRHS).*Q12;
M3p = (T.^3.*SRHS + 3.*T.^2.*S.*TRHS).*Q13;
M4p = (T.^4.*SRHS + 4.*T.^3.*S.*TRHS).*Q14;

D = -g./rho0.*(T1p + T2p + T3p + T4p + T5p+ S1p + S2p + S3p + S4p + S5p + M1p + M2p + M3p + M4p);
end