function [OUT] = Int_varROMS(Field,ingr,outgr)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% Interpolate Field onto desired Arakawa C-grid locations
%
% INPUTS:
%
% Field = input field.
%
% ingr = [hgrd vgrd] -> input grid locations.
%
% hgrd = 1 (rho), 2 (u), 3 (v), 4 (psi)
%
% vgrd = 1 (rho), 2 (w)
%
% outgrd = [hgrd vgrd] -> output grid locations.
%
% Int marks points that should be extrapolated with 0's. 
%
% For now vertical interpolation is just carried about using
% central differencing assuming a uniform vertical grid, this needs
% to be changed!!!!!!!!!!
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

[xLt,yLt,zLt,tLt] = size(Field);

hgrdI = ingr(1);
vgrdI = ingr(2);
hgrdOUT = outgr(1);
vgrdOUT = outgr(2);

%Do horizontal placement:
switch hgrdI
  case 1 % in is rho ------------------------------------------------
    switch hgrdOUT
      case 1 % out is rho
        OUT = Field;
      case 2 % out is u
        OUT = (Field(2:xLt,:,:,:)+Field(1:(xLt-1),:,:,:))/2;
      case 3 % out is v
        OUT = (Field(:,2:yLt,:,:)+Field(:,1:(yLt-1),:,:))/2;
      case 4 % out is psi
        OUT = (Field(2:xLt,2:yLt,:,:)+Field(2:xLt,1:(yLt-1),:,:)+...
               Field(1:(xLt-1),2:yLt,:,:)+Field(1:(xLt-1),1:(yLt-1),:,: ...
                                              ))/4;
    end
  case 2 % in is u --------------------------------------------------
    switch hgrdOUT
      case 1 % out is rho
        OUT = NaN*zeros(xLt+1,yLt,zLt,tLt);
        OUT(2:xLt,:,:,:) = (Field(2:xLt,:,:,:)+Field(1:(xLt-1),:,:,:))/2;
      case 2 % out is u
        OUT = Field;
      case 3 % out is v
        OUT = NaN*zeros(xLt+1,yLt-1,zLt,tLt);
        OUT(2:xLt,:,:,:) = (Field(1:(xLt-1),1:(yLt-1),:,:)+Field(2:xLt,1:(yLt-1),:,:)+...
               Field(1:(xLt-1),2:yLt,:,:)+Field(1:(xLt-1),1:(yLt-1),:,: ...
                                              ))/4;
      case 4 % out is psi
        OUT = (Field(:,2:yLt,:,:)+Field(:,1:(yLt-1),:,:))/2;
    end
  case 3 % in is v --------------------------------------------------
    switch hgrdOUT
      case 1 % out is rho
        OUT = NaN*zeros(xLt,yLt+1,zLt,tLt);
        OUT(:,2:yLt,:,:) = (Field(:,2:yLt,:,:)+Field(:,1:(yLt-1),:, ...
                                                     :))/2;
      case 2 % out is u
        OUT = NaN*zeros(xLt-1,yLt+1,zLt,tLt);
        OUT(:,2:yLt,:,:) = (Field(1:(xLt-1),1:(yLt-1),:,:)+Field(2:xLt,1:(yLt-1),:,:)+...
               Field(1:(xLt-1),2:yLt,:,:)+Field(1:(xLt-1),1:(yLt-1),:,: ...
                                              ))/4;
      case 3 % out is v
        OUT = Field;
      case 4 % out is psi
        OUT = (Field(2:xLt,:,:,:)+Field(1:(xLt-1),:,:,:))/2;
    end
  case 4 % in is psi --------------------------------------------------
    switch hgrdOUT
      case 1 % out is rho
        OUT = NaN*zeros(xLt+1,yLt+1,zLt,tLt);
        OUT(2:xLt,2:yLt,:,:) = (Field(1:(xLt-1),1:(yLt-1),:,:)+Field(2:xLt,1:(yLt-1),:,:)+...
               Field(1:(xLt-1),2:yLt,:,:)+Field(1:(xLt-1),1:(yLt-1),:,: ...
                                              ))/4;
      case 2 % out is u
        OUT = NaN*zeros(xLt,yLt+1,zLt,tLt);
        OUT(:,2:yLt,:,:) = (Field(:,1:(yLt-1),:,:)+Field(:,2:yLt,:,:))/2;
      case 3 % out is v
        OUT = NaN*zeros(xLt+1,yLt,zLt,tLt);
        OUT(2:xLt,:,:,:) = (Field(1:(xLt-1),:,:,:)+Field(2:xLt,:,:,:))/2;
      case 4 % out is psi
        OUT = Field;
    end
end

%Do vertical placement:
switch vgrdI
    case 1 % in is rho -------------------------------------------------
      switch vgrdOUT
        case 1 % out is rho
        case 2 % out is w
          OUT = NaN*zeros(xLt,yLt,zLt+1,tLt);
          OUT(:,:,2:zLt,:) = (OUT(:,:,1:(zLt-1),:)+OUT(:,:,2:zLt,:))/2;
      end
      case 2 % in is w -------------------------------------------------
      switch vgrdOUT
        case 1 % out is rho
          OUT = (OUT(:,:,1:(zLt-1),:)+OUT(:,:,2:zLt,:))/2;
        case 2 % out is w
      end
end
end