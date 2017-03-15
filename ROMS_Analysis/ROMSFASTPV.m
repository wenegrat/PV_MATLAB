function  outStruct = ROMSFASTPV(path, sliceT, pm, pn, sst, Qo, EP, SW, rho0, Cp, tx, ty,...
                                tmag, dxf, dyf, theta_s, theta_b, hc, h,zl, g, f, zmin)
r=0.58; d1=0.35; d2=23; % Solar Absorption Jerlov I
dx = squeeze(dxf(:,:,1)); dy = squeeze(dyf(:,:,1));
    % 2) Calculate z coordinates at this timestep.
    Eta = GetVarROMS(path, 0, {'zeta', '(1)'}, sliceT);

    hkpp = GetVarROMS(path, 0, {'hbls', '(1)'}, sliceT);
    [nx ny] = size(hkpp);
    
    z = compZ(path, 0, Eta, theta_s, theta_b, hc, h);
    z = z(:,:,zl);
    zwt = compZ(path, 1, Eta,  theta_s, theta_b, hc, h);
    zw(:,:,:) = zwt(:,:,zl(1):zl(end)+1);  % Save zw for later dz calc
%         zw = compZ(path, 1);
%         zw = zw(xl, yl, zl(1):zl(end)+1);
    [~, ~, nz] = size(z);
    zm = squeeze(nanmean(nanmean(z)));
    % Omega is the vertical velocity in s-coordinates
    omega = GetVarROMS(path, 0, {'omega', '(1)'}, sliceT);
    Zx = DrvROMS(pm, z, 'x'); % XXX-Is this the right type of derivative to take?
    Zy = DrvROMS(pn, z, 'y');
%     if (i>1)
%         zwt = squeeze(zw(:,:,:,i) - zw(:,:,:,i-1))./ts; % XXX- staggered back 1 timestep...
%         zwt = (zwt(:,:,2:end) + zwt(:,:,1:end-1)).*0.5;
%         zwt = 0;
%     else
        zwt = 0; % XXX-This is formally incorrect...
%     end
    
    % 3) Load Zonal Momentum Stuff
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
%     Tf(:,:,:,i)=T; %% Might be worth keeping this later for plotting purposes...
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
    % - XXX Should be adiabatically leveled first....
    Bz(:,:,:) = DrvROMS(z, B, 'z');
%     Bz(:,:,:,i) = bvf;

    
    W = zwt + U.*Zx + V.*Zy + omega;% see: https://www.myroms.org/forum/viewtopic.php?f=19&t=2139
    Wy = DrvS(pn, z, W, 'y');
    Wx = DrvS(pm, z, W, 'x');
    
%     Uf(:,:,:,i) = U; Vf(:,:,:,i) = V; Wf(:,:,:,i) = W;
    OMEGAX(:,:,:) =  - Vz;
    OMEGAY(:,:,:) =  + Uz;
    OMEGAZ(:,:,:) = repmat(f, [1 1 nz]) + Vx - Uy;
    
    Q(:,:,:) = OMEGAX(:,:,:).*Bx(:,:,:) + OMEGAY(:,:,:).*By(:,:,:)+OMEGAZ(:,:,:).*Bz(:,:,:);
    
      JAz(:,:,:) = W.*Q(:,:,:); JAx(:,:,:) = U.*Q(:,:,:); JAy(:,:,:) = V.*Q(:,:,:);

      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % MAKE THEORY SCALINGS
       % Calc Buoyancy Fluxes
      del = .1;
        alpha = squeeze(1./rho(:,:,end)).*(rho_eos(squeeze(T(:,:,end,:))+del, squeeze(S(:,:,end,:)), 0)-rho_eos(squeeze(T(:,:,end,:))-del, squeeze(S(:,:,end,:)), 0))./(2*del);
        beta =squeeze(1./rho(:,:,end)).*(rho_eos(squeeze(T(:,:,end,:)), squeeze(S(:,:,end,:))+del, 0)-rho_eos(squeeze(T(:,:,end,:)), squeeze(S(:,:,end,:))-del, 0))./(2*del);
        ssta =  squeeze(T(:,:,end))-sst(:,:);
        Qeff = Qo(:,:) - 30*ssta;
        Bo = g*alpha.*(Qeff)./(rho0*Cp) +  g*beta.*EP(:,:).*squeeze(S(:,:,end));
        BLW = g*alpha.*(Qeff-SW(:,:))./(rho0*Cp) +  g*beta.*EP(:,:).*squeeze(S(:,:,end));
        BSW = g*alpha.*SW(:,:)./(rho0*Cp);
        
        %Calc b gradients 
        ms = double(zw(:,:,1:end-1) > -repmat(hkpp, [1 1 nz]));
% pause();
        Bzh = Bz.*ms;
        bzh = NaN(nx,ny);
%         for xi=1:nx %% THIS IS PROBABLY A SLOW WAY TO DO THIS.
%             for yi=1:ny
%             ind = find(isfinite(squeeze(Bzh(xi,yi,:))), 1, 'first');
%             if ~isempty(ind)
%                 if ind<50
%                     bzh(xi,yi) = squeeze(Bz(xi,yi,ind+1));
%                 end
%             end
%             end
%         end
        
        B = findfirst(ms, 3); %% XX-Should check that this is working correctly.
        B = B+1;
        B(B>nz) = nz;
        for xi=1:nx;
            for yi =1:ny;
                if B(xi, yi)>1
                bzh(xi,yi) = Bz(xi,yi, B(xi,yi));
                end
            end
        end
%         
        ms(ms ==0) = NaN; % This has to be after the bzh calculation.

        M2 = squeeze(squeeze(nanmean(ms.*Bx(:,:,:), 3)).^2 + squeeze(nanmean(ms.*By(:,:,:), 3)).^2); %
        m = ~isfinite(M2);  
        Ma = squeeze(Bx(:,:,end-1)).^2 + squeeze(By(:,:,end-1)).^2;
        M2(m) = Ma(m); % Use uppermost gradient value if H is too small...

        % Generate scalings
        hkpp(hkpp<zmin) = zmin;

        JFT = -0.2*M2.*hkpp;

%         JBTs = repmat(f, [1 1 nt]).*Bo./hkpp;

% A more physically based way accounting for SW


        Bsw1 = r.*BSW.*(1-exp(-hkpp./d1));
        Bsw2 = (1-r).*BSW.*(1-exp(-hkpp./d2));
        Bup = BLW+Bsw1+Bsw2;
% JBTsStrat = repmat(f, [1 1 nt]).*((Bup)./Hm);
% JBTs(Bo<0) = JBTsStrat(Bo<0);
        JBTs = f.*Bup./hkpp;
        
        JENT = f.*5e-3.*bzh./(hkpp);
        JBTe = 0.15.*M2.*hkpp;
        JBT = JBTs + JBTe;
        
        % Wind
        JFW = squeeze(-tx(:,:)./(rho0*hkpp).*squeeze(By(:,:,end-2)) + ty(:,:)./(rho0*hkpp).*squeeze(Bx(:,:,end-2)));
        JFWe = 0.7*sqrt(squeeze(tmag(:,:))./rho0).^3.*f./(hkpp.^2);
        
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate Integrated Quantities
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     mask = 0.*mask; % Reset mask
    mask = (rho>1025.9) & (rho<1026.3);
%     mask = (T>18) & (T<20);
%     mask = ones(size(rho));
    dz = diff(zw, 1, 3); % Could move this out of the loop if zw is ~ constant.
    gridvol = abs(dxf.*dyf).*abs(dz);
    volt = squeeze(sum(sum(sum(mask.*gridvol))));
    outStruct.vol = nanmean(volt);
    
 
    Q(~isfinite(Q)) = 0;

    Q = Q.*mask;
    outStruct.Qat = squeeze(nansum(nansum(nansum(Q.*gridvol)))); 
    
    % Calculate Advective J Terms

    % Vert Terms
        dJAzAb = squeeze(nansum(nansum(squeeze(JAz(:,:,2).*mask(:,:,2)).*dx.*dy)));
        dJAzAt = squeeze(nansum(nansum(squeeze(JAz(:,:,end-1).*mask(:,:,end-1)).*dx.*dy)));
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
        outStruct.dJENTa = squeeze(nansum(nansum(JENT.*squeeze(mask(:,:,end-1)).*dx.*dy)));
      
        
    % Other useful thing sto keep
    outStruct.hm = squeeze(nanmean(nanmean(hkpp.*squeeze(mask(:,:,end-1)))));
    outStruct.Qm = squeeze(nanmean(nanmean(squeeze(Qo(:,:)).*squeeze(mask(:,:,end-1)))));
    outStruct.omegazs(:,:) = squeeze(OMEGAZ(:,:,end-2,:));
    outStruct.mask = mask;
    outStruct.T = T;
    outStruct.Q = Q.*gridvol;
    
%     outStruct.Bmag =sqrt(Bx.^2 + By.^2);
%     outStruct.Bo = squeeze(nanmean(nanmean(Bo.*squeeze(mask(:,:,end-1)))));
end