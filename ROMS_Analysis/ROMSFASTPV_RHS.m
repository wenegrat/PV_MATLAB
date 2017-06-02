function  [outStruct output] = ROMSFASTPV_RHS(path,pathpv, pathf,pathb, sliceT, pm, pn, sst,sss, Qo, EP, SW, rho0, Cp, tx, ty,...
                                tmag, dxf, dyf, theta_s, theta_b, hc, h,zl, g, f, zmin, ntp, ts, oldrunflag)
%%%%%%%%%%%%%%%%
% ROMSFASTPV Jacob Wenegrat (jwenegrat@stanford.edu)
%
% Function is made to be called in a parfor parallel loop
%
% INPUTS:
% path -- path to netcdf file
% pathf -- path to next netcdf file (needed for calculating dhdt)
% pathb -- path to prior netcdf file (needed for calculating dhdt)
% sliceT -- 4 cell variable denoting which slice and timestep to load.
% pm, pn -- grid metrics
% sst, sss -- climatological sst and sss for calculating relative fluxes
% Qo -- surface heat flux (interpolated to model time)
% EP -- Evap - Precip forcing (interpolated to model time)
% SW -- Shortwave radiation (for calculating penetrative component)
% ....
%%%%%%%%%%%%%%%%%

r=0.58; d1=0.35; d2=23; % Solar Absorption Jerlov I
dx = squeeze(dxf(:,:,end-1)); dy = squeeze(dyf(:,:,end-1));

    % Here because different runs have different netcdf naming conventions
    if oldrunflag
        hkppstring = 'hbls';
    else
        hkppstring = 'hbl';
    end
    
    % 2) Calculate z coordinates at this timestep.
    Eta = GetVarROMS(path, 0, {'zeta', '(1)'}, sliceT);

    hkpp = GetVarROMS(path, 0, {hkppstring, '(1)'}, sliceT); %Old runs are hbls
    
%     if (sliceT{4}(1) == ntp)
%         hkppf = GetVarROMS(pathf, 0, { hkppstring, '(1)'}, {sliceT{1}, sliceT{2}, sliceT{3}, [1 1]});
%         hkppb = GetVarROMS(path, 0, { hkppstring, '(1)'}, {sliceT{1}, sliceT{2}, sliceT{3}, [sliceT{4}(1)-1 sliceT{4}(2)-1]});
%     elseif sliceT{4}(1) == 1
%         hkppb = GetVarROMS(pathb, 0, { hkppstring, '(1)'}, {sliceT{1}, sliceT{2}, sliceT{3}, [ntp ntp]});
%         hkppf = GetVarROMS(path, 0, { hkppstring, '(1)'}, {sliceT{1}, sliceT{2}, sliceT{3}, [2 2]});
%     else % middle of one file
%         hkppb = GetVarROMS(path, 0, { hkppstring, '(1)'}, {sliceT{1}, sliceT{2}, sliceT{3}, [sliceT{4}(1)-1 sliceT{4}(2)-1]});
%         hkppf = GetVarROMS(path, 0, { hkppstring, '(1)'}, {sliceT{1}, sliceT{2}, sliceT{3}, [sliceT{4}(1)+1 sliceT{4}(2)+1]});
%     end
%     dhdt = (hkppf-hkppb)./(2*ts);
%     dhdt(dhdt<0) = 0;
    
    [nx ny] = size(hkpp);
    
    z = compZ(path, 0, Eta, theta_s, theta_b, hc, h);
    z = z(:,:,zl);
    zwt = compZ(path, 1, Eta,  theta_s, theta_b, hc, h);
    zw(:,:,:) = zwt(:,:,zl(1):zl(end)+1);  % Save zw for later dz calc
    [~, ~, nz] = size(z);
    zm = squeeze(nanmean(nanmean(z)));
    % Omega is the vertical velocity in s-coordinates
%     omega = GetVarROMS(path, 0, {'omega', '(1)'}, sliceT);
    Zx = DrvROMS(pm, z, 'x'); % XXX-Is this the right type of derivative to take?
    Zy = DrvROMS(pn, z, 'y');
%     if (i>1)
%         zwt = squeeze(zw(:,:,:,i) - zw(:,:,:,i-1))./ts; % XXX- staggered back 1 timestep...
%         zwt = (zwt(:,:,2:end) + zwt(:,:,1:end-1)).*0.5;
%         zwt = 0;
%     else
        zwt = 0; % XXX-This is formally incorrect...
%     end
    
    % 3) Load Zonal Vels
    U = GetVarROMS(path, 0, {'u', '(1)'}, sliceT, 1);
    Uu = GetVarROMS(path, 0, {'u', '(1)'}, sliceT, 0);
    % Advective Terms
    Uz = DrvROMS(z, U, 'z');
    Uy = DrvS(pn,z, Uu, 'y', 2);
    
    % 4) Load Meridional Velocity
    V = GetVarROMS(path, 0, {'v', '(1)'}, sliceT);
    Vv = GetVarROMS(path, 0, {'v', '(1)'}, sliceT,0);

    % Advective V
    Vz = DrvROMS(z, V, 'z');
    Vx = DrvS(pm,z, Vv, 'x', 3);
   
    
    % 5) Load Buoyancy Terms
    T = GetVarROMS(path, 0, {'temp', '(1)'}, sliceT);
    S = GetVarROMS(path, 0, {'salt', '(1)'}, sliceT);
    rho = rho_eos(T, S, 0); % CROCO function (checked for consistency 1/12/17)
    B = -g*rho./rho0; % In-Situ B
    % This calculates adiabatically referenced N^2    
    ztemp = permute(z, [3 1 2]);
    zwtemp = permute(squeeze(zw(:,:,1:end-1)), [3 1 2]);
    Ttemp = permute(T, [3 1 2]);
    Stemp = permute(S, [3 1 2]);
    [rhoa, bvf] = rho_eos(Ttemp,Stemp,ztemp, zwtemp, g, rho0);
    bvf = permute(bvf, [2 3 1]); % N^2 
    

    Bx(:,:,:) = DrvS(pm, z, B, 'x');
    By(:,:,:) = DrvS(pn, z, B, 'y');
    Bz(:,:,:) = DrvROMS(z, B, 'z'); % Note, could use either definition...not sure if one is to be preferred.
%     Bz(:,:,:) = bvf; disp('Alternte N2 def');

    disp('Warning: No W loaded'); % The high res runs don't have W output...
%     W = zwt + U.*Zx + V.*Zy + omega;% see: https://www.myroms.org/forum/viewtopic.php?f=19&t=2139
%     Wy = DrvS(pn, z, W, 'y');
%     Wx = DrvS(pm, z, W, 'x');
    
%     Uf(:,:,:,i) = U; Vf(:,:,:,i) = V; Wf(:,:,:,i) = W;
    OMEGAX(:,:,:) =  - Vz;
    OMEGAY(:,:,:) =  + Uz;
    OMEGAZ(:,:,:) = repmat(f, [1 1 nz]) + Vx - Uy;
    
    % Define PV
    Q(:,:,:) = OMEGAX(:,:,:).*Bx(:,:,:) + OMEGAY(:,:,:).*By(:,:,:)+OMEGAZ(:,:,:).*Bz(:,:,:);
    
    % Advective fluxes
    JAz(:,:,:) = zeros(size(Q)); %W.*Q(:,:,:); 
    JAx(:,:,:) = U.*Q(:,:,:); JAy(:,:,:) = V.*Q(:,:,:);

    dz = diff(zw, 1, 3); % Could move this out of the loop if zw is ~ constant.

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % MAKE THEORY SCALINGS
       % Calc Buoyancy Fluxes
        del = .1;
        alphaf = squeeze(1./rho(:,:,:)).*(rho_eos(squeeze(T(:,:,:,:))+del, squeeze(S(:,:,:,:)), 0)-rho_eos(squeeze(T(:,:,:,:))-del, squeeze(S(:,:,:,:)), 0))./(2*del);
        alpha = squeeze(alphaf(:,:,end));
        betaf = squeeze(1./rho(:,:,:)).*(rho_eos(squeeze(T(:,:,:,:)), squeeze(S(:,:,:,:))+del, 0)-rho_eos(squeeze(T(:,:,:,:)), squeeze(S(:,:,:,:))-del, 0))./(2*del);
        beta =  squeeze(betaf(:,:,end));
        ssta =  squeeze(T(:,:,end))-sst(:,:);
        Qeff = Qo(:,:) - 30*ssta;
        Teff = Qeff./(rho0.*Cp);
        Seff = EP.*squeeze(S(:,:,end)) - 1./(90.*86400).*(squeeze(S(:,:,end)) - sss).*squeeze(dz(:,:,end));
        
        Bo = g*alpha.*Teff +  g*beta.*Seff;
        BLW = g*alpha.*(Teff-SW(:,:)./(rho0*Cp)) +  g*beta.*Seff;
        BSW = g*alpha.*SW(:,:)./(rho0*Cp);
        
        %Calc b gradients 
        ms = double(zw(:,:,1:end-1) >= -repmat(hkpp, [1 1 nz]));
        bzh = NaN(nx,ny);
        akh  = bzh;
        bdiff = bzh;
        Akt = GetVarROMS(path, 0, {'AKt', '(1)'}, sliceT);
       
        BI = findfirst(ms, 3); %%Need to confirm this function is working correctly.
        BI(BI>nz) = nz;
        for xi=1:nx;
            for yi =1:ny;
                if BI(xi, yi)>1
                bzh(xi,yi) = Bz(xi,yi, BI(xi,yi));
                akh(xi,yi) = Akt(xi,yi,BI(xi,yi));
                bmean = nanmean(B(xi,yi,BI(xi,yi):end));
                bdiff(xi,yi) = bmean - B(xi,yi, BI(xi, yi)-1); % XXX-CHeck Sensitivity to how this is calculated.
                end
            end
        end
%         
        ms(ms ==0) = NaN; % This has to be after the bzh calculation.
        Bya = squeeze(nanmean(ms.*By(:,:,:), 3)); %Average over the boundary layer
        Bxa = squeeze(nanmean(ms.*Bx(:,:,:), 3));
        
        M2 = Bxa.^2 + Bya.^2; %
        m = ~isfinite(M2); 
        Ma = squeeze(Bx(:,:,end-1)).^2 + squeeze(By(:,:,end-1)).^2;
        M2(m) = Ma(m); % Use uppermost gradient value if H is too small...
        N2 = squeeze(nanmean(Bz.*ms,3));

        % Generate scalings
        hkpp(hkpp<zmin) = zmin;
        JFT = -0.2*M2.*hkpp;


% A more physically based way accounting for SW


        Bsw1 = r.*BSW.*(1-exp(-hkpp./d1));
        Bsw2 = (1-r).*BSW.*(1-exp(-hkpp./d2));
        Bup = BLW+Bsw1+Bsw2;
        Bup = Bo; %Not accounting for penetrating radiation.

        
        JBTs = f.*Bup./hkpp;
        
%         JENT = 0.*f.*akh.*bzh./(hkpp) + f.*dhdt.*bdiff./hkpp; %f.*dhdt.*bzh;
        JBTe = 0.15.*M2.*hkpp;
        JBT = JBTs + JBTe;
        
        % Wind
        dele = 0.4.*sqrt(abs(tx+1i*ty)./rho0)./f;
        dele = hkpp;
        JFW = squeeze(-tx(:,:)./(rho0*dele).*Bya + ty(:,:)./(rho0*dele).*Bxa);
        JFWe = 0.7*sqrt(squeeze(tmag(:,:))./rho0).^3.*f./(hkpp.^2);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Integrated Quantities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     mask = 0.*mask; % Reset mask
    mask = (rho>1025.9) & (rho<1026.3);
    mask = (rho>1024.9) & (rho<1026.9); disp('Alternate Mask');
%     mask = (T>18) & (T<20);
%     mask = ones(size(rho)); disp('Using full volume mask!');
    gridvol = abs(dxf.*dyf).*abs(dz);
    volt = squeeze(sum(sum(sum(mask.*gridvol))));
    outStruct.vol = nanmean(volt);
    
 
    Q(~isfinite(Q)) = NaN;

    Q = Q.*mask;
    outStruct.Qat = squeeze(nansum(nansum(nansum(Q.*gridvol.*rho))))./1027; 
    
    % Calculate Advective J Terms

    % Vert Terms
        dJAzAb = squeeze(nansum(nansum(squeeze(JAz(:,:,2).*mask(:,:,2)).*dx.*dy)));
%         dJAzAt = squeeze(nansum(nansum(squeeze(JAz(:,:,end-1).*mask(:,:,end-1)).*dx.*dy)));
        outStruct.dJAzA =  - dJAzAb; % NO ADV FLUX OUT SURFACE...
    % Zonal Terms
        dJAx_l = squeeze(nansum(nansum(squeeze(JAx(2,:,:).*mask(2,:,:)).*squeeze(dyf(2,:,:).*dz(2,:,:)))));
        dJAx_r = squeeze(nansum(nansum(squeeze(JAx(end-1,:,:).*mask(end-1,:,:)).*squeeze(dyf(end-1,:,:).*dz(end-1,:,:)))));
        outStruct.dJAxA = dJAx_r - dJAx_l;
    % Meridional Terms
        dJAy_f = squeeze(nansum(nansum(squeeze(JAy(:,2,:).*mask(:,2,:)).*squeeze(dxf(:,2,:).*dz(:,2,:)))));
        dJAy_b = squeeze(nansum(nansum(squeeze(JAy(:,end-1,:).*mask(:,end-1,:)).*squeeze(dxf(:,end-1,:).*dz(:,end-1,:)))));
        outStruct.dJAyA = dJAy_b - dJAy_f;
        
   % Calculate Scaling Terms
        outStruct.dJBTA = squeeze(nansum(nansum(JBT.*squeeze(mask(:,:,end-1)).*dx.*dy)));
        outStruct.dJFTA = squeeze(nansum(nansum(JFT.*squeeze(mask(:,:,end-1)).*dx.*dy)));        
        outStruct.dJFWA = squeeze(nansum(nansum(JFW.*squeeze(mask(:,:,end-1)).*dx.*dy)));
        outStruct.dJBTAs = squeeze(nansum(nansum(JBTs.*squeeze(mask(:,:,end-1)).*dx.*dy)));
        outStruct.dJBTAe = squeeze(nansum(nansum(JBTe.*squeeze(mask(:,:,end-1)).*dx.*dy)));
        outStruct.dJBFWe = squeeze(nansum(nansum(JFWe.*squeeze(mask(:,:,end-1)).*dx.*dy)));
%         outStruct.dJENTa = squeeze(nansum(nansum(JENT.*squeeze(mask(:,:,end-1)).*dx.*dy)));
      
        
    % Other useful thing sto keep
    outStruct.hm = squeeze(nanmean(nanmean(hkpp.*squeeze(mask(:,:,end-1)))));
    outStruct.h = squeeze(hkpp);
    outStruct.Bo = Bo;
%     outStruct.STRAIN = squeeze(Uy+Vx);
%     outStruct.M2 = M2;
%     outStruct.N2 = N2;
    outStruct.Qm = squeeze(nanmean(nanmean(squeeze(Qo(:,:)).*squeeze(mask(:,:,end-1)))));
    outStruct.omegazs(:,:) = squeeze(OMEGAZ(:,:,end-2,:));
%     outStruct.OMEGAZ = OMEGAZ;
%     outStruct.OMEGAX = OMEGAX;
%     outStruct.OMEGAY = OMEGAY;
    outStruct.JFT = JFT;
    outStruct.JBTs = JBTs;
    outStruct.JBTe = JBTe;
    outStruct.JFW = JFW;
%     outStruct.JENT = JENT;
    outStruct.Bx = Bx;
    outStruct.By = By;
    outStruct.Bz = Bz;
%     outStruct.alpha = alphaf;
%     outStruct.beta = betaf;
    outStruct.mask = mask;
%     outStruct.T = T;
%     outStruct.Q = Q;%.*gridvol;
    outStruct.dz = dz;
%     outStruct.Rho = rho;
    outStruct.dz = abs(dz);
    outStruct.negQ = squeeze(nansum(nansum(nansum(Q<0))));
    
    
    %% RHS TERMS
    % U RHS
URHS = GetVarROMS(pathpv, 0, {'u_rhs', '(1)'}, sliceT);
% V RHS
VRHS = GetVarROMS(pathpv, 0, {'v_rhs', '(1)'}, sliceT);

% T RHS
TRHS = GetVarROMS(pathpv, 0, {'temp_rhs', '(1)'}, sliceT);

% S RHS
SRHS = GetVarROMS(pathpv, 0, {'salt_rhs', '(1)'}, sliceT);



JFZ = Bx.*VRHS - By.*URHS;
JFX = -Bz.*VRHS;
JFY = Bz.*URHS;

% D   = g*alpha.*(TRHS)   + g*beta.*(SRHS);
% D = -g./1027.*rho_eos(TRHS, SRHS, 0);
D = bRHS(TRHS, SRHS, T, S, 1027.4, g);

% magb = sqrt(Bx.^2 + By.^2 + Bz.^2);
% FACT = (Q.*Bz./(magb.^2)); % See HM90 Eq. 4.1 and 4.2
JDZ = -(OMEGAZ).*D;
JDX = -OMEGAX.*D;
JDY = -OMEGAY.*D;

% XPROJ = Bx./magb; 
% YPROJ = By./magb;
% ZPROJ = Bz./magb;

output.JFZ = squeeze(JFZ(:,:,end-1));
output.JDZ = squeeze(JDZ(:,:,end-1));

output.Jfa = squeeze(nansum(nansum(squeeze(JFZ(:,:,end-1).*mask(:,:,end-1).*dxf(:,:,end-1).*dyf(:,:,end-1)))));
output.Jda = squeeze(nansum(nansum(squeeze(JDZ(:,:,end-1).*mask(:,:,end-1).*dxf(:,:,end-1).*dyf(:,:,end-1)))));

[nx ny nz] = size(D);
% X-Fluxes
xl = 2; xr = nx-1;
jfxl = squeeze(nansum(nansum(squeeze(JFX(xl,:,:).*mask(xl,:,:).*dyf(xl,:,:).*dz(xl,:,:)))));
jfxr = squeeze(nansum(nansum(squeeze(JFX(xr,:,:).*mask(xr,:,:).*dyf(xr,:,:).*dz(xr,:,:)))));
jdxl = squeeze(nansum(nansum(squeeze(JDX(xl,:,:).*mask(xl,:,:).*dyf(xl,:,:).*dz(xl,:,:)))));
jdxr = squeeze(nansum(nansum(squeeze(JDX(xr,:,:).*mask(xr,:,:).*dyf(xr,:,:).*dz(xr,:,:)))));

% Y-Fluxes
yf = 2; yb = ny-1;
jfyf = squeeze(nansum(nansum(squeeze(JFY(:,yf,:).*mask(:,yf,:).*dxf(:,yf,:).*dz(:,yf,:)))));
jfyb = squeeze(nansum(nansum(squeeze(JFY(:,yb,:).*mask(:,yb,:).*dxf(:,yb,:).*dz(:,yb,:)))));
jdyf = squeeze(nansum(nansum(squeeze(JDY(:,yf,:).*mask(:,yf,:).*dxf(:,yf,:).*dz(:,yf,:)))));
jdyb = squeeze(nansum(nansum(squeeze(JDY(:,yb,:).*mask(:,yb,:).*dxf(:,yb,:).*dz(:,yb,:)))));

 output.Jfa = output.Jfa + (jfxr - jfxl)+(jfyb-jfyf);
 output.Jda = output.Jda + (jdxr - jdxl)+(jdyb-jdyf);
end