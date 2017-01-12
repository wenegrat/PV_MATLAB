function [OUT] = DrvROMS(metric,Field,TYPE,varargin)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function takes the central difference derivative of Field on
% the grid grd in the direction direc and then interpolates the
% result back onto the rho grid.
%
% INPUTS:
% TYPE = 'x' - x (always 1st dimension derivative)
%      = 'y' - y (always 2nd dimension derivative)
%      = 'z' - z (always 3rd dimension derivative)
%
% metric if TYPE = 'x' pn (2D horizontal)
% metric if TYPE = 'y' pm (2D horizontal)
% metric if TYPE = 'z' z_rho or z_w (3D)
%
% Field = the field to take the derivative on.
%
% varargin =  or ,hgrd or hgrd,vgrd
%
% hgrd = 1 (rho -default), 2 (u), 3 (v), 4 (psi).
%
% vgrd = 1 (rho -default), 2 (w).
%
%---------------------------------------------------------------------
%
% Dependencies; Int.
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------

% Get grid types:
    hgrdD = 1;
    vgrdD = 1;
    if (nargin == 4)
        hgrdD = varargin{1};
    elseif (nargin == 5)
        hgrdD = varargin{1};
        vgrdD = varargin{2};
    end
    [d1,d2,d3,d4] = size(Field);
    
%     metric = Int_varROMS(metric,[1 1],[hgrdD 1]);
% Take the derivative:
    switch TYPE
      case 'x' % x derivative metric = pn -----------------------------------
        if (d1 < 3)
            warning(['Can''t take x derivative for less than 3 x ' ...
                     'points!']);
        elseif (d2 < 2 && hgrdD == 3)
            warning(['Can''t take x derivative on v grid for less ' ...
                     'than 2 u points']);
        end
        
        switch hgrdD
          case 1 % rho - interp pn.
            dist = repmat(1./metric(2:end,:)+1./metric(1:(end-1),:),[1 1 d3 d4])/2;            
            OUT = (Field(2:d1,:,:,:)-Field(1:(d1-1),:,:,:))./dist;
            hgrdD = 2;
          case 2 % u grid - no change to pn
            OUT = NaN*zeros(d1+1,d2,d3,d4);
            dist = repmat(1./metric,[1 1 d3 d4]);
            OUT(2:d1,:,:,:) = (Field(2:d1,:,:,:)-Field(1:(d1-1),:,: ...
                                                       ,:))./dist(2:d1,:,:,:);
            hgrdD = 1;
          case 3 % v grid - interp pn.
            dist = repmat(1./metric(2:end,2:end)+1./metric(2:end,1:(end-1))+...
                          1./metric(1:(end-1),2:end)+1./metric(1:(end-1),1:(end-1))...
                          ,[1 1 d3 d4])/4;
            OUT = (Field(2:d1,:,:,:)-Field(1:(d1-1),:,:,:))./dist;
            hgrdD = 4;
          case 4 % psi grid - interp pn.
            OUT = NaN*zeros(d1+1,d2,d3,d4);
            dist = repmat(1./metric(:,2:end)+1./metric(:,1:(end-1)),[1 1 d3 d4])/2;
            OUT(2:d1,:,:,:) = (Field(2:d1,:,:,:)-Field(1:(d1-1),:,:,:))./dist(2:d1,:,:,:);
            hgrdD = 3;
        end
      case 'y' % y derivative metric = pm ----------------------------------
        if (d2 < 3)
            warning(['Can''t take y derivative for less than 3 y ' ...
                     'points!']);
        elseif (d1 < 1 && hgrdD == 2)
            warning(['Can''t take y derivative on u grid for less ' ...
                     'than 2 u points']);
        end
        
        switch hgrdD
          case 1 % rho - interp pm.
            dist = repmat(1./metric(:,2:end)+1./metric(:,1:(end-1)),[1 1 d3 d4])/2;
            OUT = (Field(:,2:d2,:,:)-Field(:,1:(d2-1),:,:))./dist;
            hgrdD = 3;
          case 2 % u grid - no change to pm.
            dist = repmat(1./metric(2:end,2:end)+1./metric(2:end,1:(end-1))+...
                          1./metric(1:(end-1),2:end)+1./metric(1:(end-1),1:(end-1))...
                          ,[1 1 d3 d4])/4;
            OUT = (Field(:,2:d2,:,:)-Field(:,1:(d2-1),:,:))./dist;
            hgrdD = 4;
          case 3 % v grid - interp pm.
            OUT = NaN*zeros(d1,d2+1,d3,d4);
            dist = repmat(1./metric,[1 1 d3 d4]);            
            OUT(:,2:d2,:,:) = (Field(:,2:d2,:,:)-Field(:,1:(d2-1),:,:))./dist(:,2:d2,:,:);
            hgrdD = 1;
          case 4 % psi grid - interpolate pm.
            dist = repmat(1./metric(2:end,:)+1./metric(1:(end-1),:),[1 1 d3 d4])/2;            
            OUT = NaN*zeros(d1,d2+1,d3,d4);
            OUT(:,2:d2,:,:) = (Field(:,2:d2,:,:)-Field(:,1:(d2-1),:,:))./dist(:,2:d2,:,:);
            hgrdD = 2;
        end
      case 'z' % z derivative
               % ---------------------------------------------
        if (d3 < 3)
            warning(['Can''t take z derivative for less than 3 z ' ...
                     'points!']);
        end
        
        %adjust vertical metric for horizontal grid:
        metric = Int_varROMS(metric,[1 1],[hgrdD 1]);
        
        switch vgrdD
          case 1 % rho
            OUT = NaN*zeros(d1,d2,d3+1,d4);
            OUT(:,:,2:d3,:) = (Field(:,:,2:d3,:)-Field(:,:,1:(d3-1),:))./ ...
                  (metric(:,:,2:d3,:)-metric(:,:,1:(d3-1),:));
          case 2 % w
            OUT = (Field(:,:,2:d3,:)-Field(:,:,1:(d3-1),:))./ ...
                  (metric(:,:,2:d3,:)-metric(:,:,1:(d3-1),:));
        end
        vgrdD = 3-vgrdD;
    end
    % Put on rho grid:
    OUT = Int_varROMS(OUT,[hgrdD vgrdD],[1 1]);
end