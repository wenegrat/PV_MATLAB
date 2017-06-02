function [var] = GetVarROMS(hname,dname,VarOp,slice, varargin);
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function returns a specified slice (slice) of the specified
% variable (VarOp) on the rho grid from the ROMS generated files hname
% and dname.
%
% INPUTS:
%
% hname = _his.nc history file generated by ROMS.
%
% dname = _dia.nc diagnostics file generated by ROMS.
%
% VarOp = variable definition string cell.
% has form:
% VarOp = {'Var1'} - only 1 variable
% VarOp = {'Var1' 'Var2'} - defaults to addition. 
% VarOp = {'Var1' 'Var2' 'Operation'}
%
% slice = cell defining the slice to take, for this it is always in
% the form;
%
% {[x1 x2],[y1 y2],[d1 d2],[t1 t2]}; slice using range x1:x2, y1:y2,
% d1:d2, t1:t2. if any of [x1 x2] = 0 then take entire range.
%
% {[x1 x2],[y1 y2],[t1 t2]}; 2D slice for a 3D variable (no depth).
%
% {[x1 x2],[y1 y2],rho,[t1 t2]}; 2D isopycnal slice; (rho positive
% real)
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%--------------Get file information-----------------------------------
% acstring = 'Single';

if (nargin == 5)
        intFlag = varargin{1};
else
    intFlag = 1;
end

ncid = netcdf.open(hname,'NC_NOWRITE');
if (strfind(hname,'avg'))
    avgfile = logical(1);
else
    avgfile = logical(0);
end
if (dname ~= 0)
    ncidD = netcdf.open(dname,'NC_NOWRITE');
end
%extras file name:
ename = [hname(1:(end-3)) '_ext.nc'];
if (exist(ename) == 2)
    ncidE = netcdf.open(ename,'NC_NOWRITE');
end

%Lengths of file:
[tmp,xLhis] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid, 'xi_rho'));
[tmp,yLhis] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid, 'eta_rho'));
[tmp,zLhis] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid, 's_rho'));
str = 'ocean_time';
try
[tmp,tLhis] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid,str));
catch exception
    if (strcmp(exception.identifier,['MATLAB:imagesci:netcdf:' ...
                            'libraryFailure']))
        str = 'clm_time';
try
[tmp,tLhis] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid,str));
catch exception
    if (strcmp(exception.identifier,['MATLAB:imagesci:netcdf:' ...
                            'libraryFailure']))
        str = 'time'; %Was temp_time, JOW 1/10/16
[tmp,tLhis] = netcdf.inqDim(ncid,netcdf.inqDimID(ncid,str));
    end
end
    end
end

%--------------Get domain---------------------------------------------
%determine type of slice:
if (length(slice) == 4)
    if (length(slice{3}) == 1)
        if (slice{3} == 0)
            vtype = 1;
        else
            vtype = 4; %isopycnal slice
        end
    else
    vtype = 1; % normal slice of 4D variable
    end
elseif (length(slice) == 3)
    vtype = 2; % normal slice of 3D variable
% $$$ elseif (length(slice) == 2)
% $$$     vtype = 3; % normal slice of 2D variable
end
if (slice{1} == 0)
    xi = 1;
    xL = xLhis;
else
xi = slice{1}(1);
xL = slice{1}(2)-xi+1;
end

if (slice{2} == 0)
    yi = 1;
    yL = yLhis;
else
yi = slice{2}(1);
yL = slice{2}(2)-yi+1;
end

if (slice{end} == 0)
    ti = 1;
    tL = tLhis;
else
    ti = slice{end}(1);
    tL = slice{end}(2)-ti+1;
end

if (vtype == 4)
    rhoISP = slice{3};
    zi = 1;
    zL = zLhis;
elseif (vtype == 1)
if (slice{3} == 0)
    zi = 1;
    zL = zLhis;
else
    zi = slice{3}(1);
    zL = slice{3}(2)-zi+1;
end
end
%------------Build operation string-----------------------------------

%If no operation string use default:
if (length(strfind(VarOp{end}, '(1)')) == 0)
    numvars = length(VarOp);
    VarOp{numvars+1} = '(1)';
    for i=2:numvars
        VarOp{numvars+1} = [VarOp{numvars+1} ' + (' num2str(i) ...
                            ')'];
    end
end   

if (length(VarOp) == 2)
if (strcmp(VarOp{1},'PV') & strcmp(VarOp{2},'(1)'))
    VarOp = {'u' 'v' 'f' 'b' ['((3)+Dx(2)-Dy(1)).*Dz(4)-Dz(2).*' ...
                        'Dx(4)+Dz(1).*Dy(4)']};
elseif (strcmp(VarOp{1},'PVv')  & strcmp(VarOp{2},'(1)'))
    VarOp = {'u' 'v' 'f' 'b' '((3)+Dx(2)-Dy(1)).*Dz(4)'};
elseif (strcmp(VarOp{1},'PVh')  & strcmp(VarOp{2},'(1)'))
    VarOp = {'u' 'v' 'f' 'b' '-Dz(2).*Dx(4)+Dz(1).*Dy(4)'};
end
end
% $$$ for i=1:length(VarOp)
% $$$     VarOp{i}
% $$$ end

GetVarStr = VarOp{end}; %Operation string.

%Custom simple variables:
% $$$ while (vari ~= length(VarOp))
% $$$     if (strcmp(VarOp{vari},'PV'));
% $$$        for j = 1:(length(VarOp)-1) 
% $$$     end
% $$$     if (strcmp(VarOp{vari},'PVh'))
% $$$     
% $$$     end
% $$$     if (strcmp(VarOp{vari},'PVv'))
% $$$         
% $$$     end
% $$$     

for i = 1:(length(VarOp)-1)
if (strcmp(VarOp{i},'b')) % Buoyancy
        VarOp{i} = 'rho';
        GetVarStr = strrep(GetVarStr,['(' num2str(i) ')'], ...
                           ['(' num2str(i) ')*(-0.0096)']);
end
end
% $$$ end

%Exxtend range if needed:
uext = 0;
vext = 0;
if (xL == 1 | yL == 1) %need to extend x or y range? 
    if (~(strcmp(VarOp{i},'z_rho') | strcmp(VarOp{i},'SF') | ...
          strcmp(VarOp{i},'TIVbot') | strcmp(VarOp{i},'TIVtop') | ...
          strcmp(VarOp{i},'z_w')))
    for i=1:(length(VarOp)-1)
    if (length(strfind(VarOp{i}, 'u_')) > 0 | ...
        length(strfind(VarOp{i}, 'v_')) > 0 | ...
        length(strfind(VarOp{i}, 'temp_')) > 0 | ...
        length(strfind(VarOp{i}, 'salt_')) > 0) %normal
        ncidC = ncidD;
    elseif (length(strfind(VarOp{i},'SF')) > 0 | ...
            length(strfind(VarOp{i},'TIVtop')) > 0 | ...
            length(strfind(VarOp{i},'TIVbot')) > 0 | ...
            length(strfind(VarOp{i},'z_rho')) > 0) %Custom
                                                   %extra
        ncidC = ncidE;
    else %Diagnostic variable. 
        ncidC = ncid;
    end

    [hgrd,vgrd,Vdims] = grd_var(ncidC,VarOp{i});
    if (hgrd == 2)
        uext = 1;
    elseif (hgrd == 3)
        vext = 1;
    end
    end
    
    if ((xL == 1 & uext == 1) | (xL == 1 & vtype == 4)) %extend domain in x
        xi = xi-1;
        xL = xL+2;
        uext = 1;
    else
        uext = 0;
    end
    if ((yL == 1 & vext == 1) | (yL == 1 & vtype == 4)) %extend domain in y
        yi = yi-1;
        yL = yL+2;
        vext = 1;
    else
        vext = 0;
    end 
    end
end
%Extend z range if needed:
zext = 0;
textn = 0;
if (vtype ~= 2)
    if (zL == 1 & length(strfind(VarOp{end},'Dz'))>0)
        zi = zi-1;
        zL = zL+2;
        zext = 1;
    end
end


%Test for horizontal derivatives:
if (length(strfind(GetVarStr,'Dx'))>0)
    Dxswitch = 1;
else
    Dxswitch = 0;
end
if (length(strfind(GetVarStr,'Dy'))>0)
    Dyswitch = 1;
else
    Dyswitch = 0;
end
if (length(strfind(GetVarStr,'Dz'))>0)
    Dzswitch = 1;
else
    Dzswitch = 0;
end
Dzrswitch = 0;
Dzwswitch = 0;    

%Get raw variable strings:
for i=1:(length(VarOp)-1)
    %Custom variables ----------------------------------------
    if (strcmp(VarOp{i},'z_rho')) %z_rho
        Dzrswitch = 1;
        GetVarStr = strrep(GetVarStr,['(' num2str(i) ')'], ...
                           'z_rho');
    elseif (strcmp(VarOp{i},'z_w')) %z_w
        Dzwswitch = 1;
        GetVarStr = strrep(GetVarStr,['(' num2str(i) ')'], ...
                           'z_w');
        
    %Normal variables ----------------------------------------
    else
%         if (length(strfind(VarOp{i}, 'u_')) > 0 | ...
%             length(strfind(VarOp{i}, 'v_')) > 0 | ...
%             length(strfind(VarOp{i}, 'temp_')) > 0 | ...
%             length(strfind(VarOp{i}, 'salt_')) > 0) %Diagnostic
% %             ncidC = ncidD;
%         elseif (length(strfind(VarOp{i},'SF')) > 0 | ...
%                 length(strfind(VarOp{i},'TIVtop')) > 0 | ...
%                 length(strfind(VarOp{i},'TIVbot')) > 0) %Custom
%                                                         %extra
%             ncidC = ncidE;
%         else
            ncidC = ncid;
%         end

        [hgrd,vgrd,Vdims] = grd_var(ncidC,VarOp{i});

        %GetVar extract:
        GetVar = ['netcdf.getVar(' num2str(ncidC) ',netcdf.inqVarID(' num2str(ncidC) ',''' ...
                  VarOp{i} '''),' Dom(ncidC,VarOp{i}) ...
                  ',''double'')'];
        if (strcmp(VarOp{i},'rho'))
            GetVar = ['(sw_dens0((' ...
                      'netcdf.getVar(' num2str(ncidC) ',netcdf.inqVarID(' ...
                      num2str(ncidC) ',''salt''),' Dom(ncidC,'salt') ...
                  ',''double'')),(' ...
                      'netcdf.getVar(' num2str(ncidC) ',netcdf.inqVarID(' ...
                      num2str(ncidC) ',''temp''),' Dom(ncidC,'temp') ...
                  ',''double'')))-1000)'];
        end                      
        if (strcmp(VarOp{i},'BRHS'))
            GetVar = ['(sw_dens0((' ...
                      'netcdf.getVar(' num2str(ncidD) ',netcdf.inqVarID(' ...
                      num2str(ncidD) ',''salt_vdiff''),' Dom(ncidD,'salt_vdiff') ...
                  ',''double'')),(' ...
                      'netcdf.getVar(' num2str(ncidD) ',netcdf.inqVarID(' ...
                      num2str(ncidD) ',''temp_vdiff''),' Dom(ncidD,'temp_vdiff') ...
                  ',''double''))))'];
        end
 
               
 
        %Interpolate for diagnostics variable:
        if (dname ~= 0 & ~avgfile)
            if (ncidC == ncidD)
                if (ti > 1 & (ti+tL-1) < tLhis) %No problems with end points.
                    GetVar = ['((' GetVar ' + ' strrep(GetVar,'ti-1','ti-2') ...
                              ')/2)']
                elseif (ti<=1)
                    warning(['Haven''t written this yet!']);
                elseif (ti+tL-1>=tLhis)
                    tL = tL-1;
                    textn = 1;
                end
            end
        end
        
        %Do repmat for 2D variables:
        if ((vtype == 4 | vtype == 1) & Vdims == 2) %2D variable in
                                                    %4D plot.
            GetVar = ['repmat(' GetVar ',[1 1 zL tL])'];
        elseif (vtype == 2 & Vdims == 2) %2D variable in 3D plot.
            GetVar = ['repmat(' GetVar ',[1 1 tL])'];
        end
        
        %Taking x derivative:
        GetVarStr = strrep(GetVarStr,['Dx(' num2str(i) ')'],...
                    ['Drv(pn,' GetVar ',''x'',' num2str(hgrd) ',' ...
                     num2str(vgrd) ')']);
        %Taking y derivative:
        GetVarStr = strrep(GetVarStr,['Dy(' num2str(i) ')'],...
                    ['Drv(pm,' GetVar ',''y'',' num2str(hgrd) ',' ...
                     num2str(vgrd) ')']);
        %Taking z derivative:
        if (Dzswitch == 1)
        if (vgrd == 1)
        GetVarStr = strrep(GetVarStr,['Dz(' num2str(i) ')'],...
                    ['Drv(z_rho,' GetVar ',''z'',' num2str(hgrd) ',' ...
                     num2str(vgrd) ')']);
        Dzrswitch = 1;
        else
        GetVarStr = strrep(GetVarStr,['Dz(' num2str(i) ')'],...
                    ['Drv(z_w,' GetVar ',''z'',' num2str(hgrd) ',' ...
                     num2str(vgrd) ')']);
        Dzwswitch = 1;
        end
        end
        %Otherwise interpolate variable:
        if (intFlag)
        GetVarStr = strrep(GetVarStr,['(' num2str(i) ')'], ['Int(' ...
                            GetVar ',[' num2str(hgrd) ' ' num2str(vgrd) ...
                            '],[1 1])']);
        else
        GetVarStr = strrep(GetVarStr,['(' num2str(i) ')'], GetVar);
        end
    end
end
    %pm, pn and z_rho are easy to deal with:
    
    if (Dxswitch == 1)
        pn = 0;
    eval(['pn = netcdf.getVar(ncid,netcdf.inqVarID(ncid,''pn''),' ...
             Dom(ncid,'pm') ',''double'');']);
    end
    if (Dyswitch == 1)
        pm = 0;
    eval(['pm = netcdf.getVar(ncid,netcdf.inqVarID(ncid,''pm''),' ...
             Dom(ncid,'pm') ',''double'');']);
    end
    if (Dzrswitch == 1 | vtype == 4)
        if (exist(ename) == 2) %use z_rho from extras file if
                               %exists. 
        eval(['z_rho = netcdf.getVar(ncidE,netcdf.inqVarID(ncidE,''z_rho''),' ...
              Dom(ncidE,'z_rho') ',''double'');']);
        else
            z_rho = zeros(xL,yL,zL,tL);
            for tind=1:tL
                ztmp = get_depths(hname,1,0,tind+ti-1,[xi yi xL ...
                                    yL]);
                z_rho(:,:,:,tind) = ztmp(:,:,zi:(zi+zL-1));
            end 
        end
    end
    if (Dzwswitch == 1)
        z_w = zeros(xL,yL,zL+1,tL);
        for tind=1:tL
            ztmp = get_depths(hname,5,0,tind+ti-1,[xi yi xL yL]);
            z_w(:,:,:,tind) = ztmp(:,:,zi:(zi+zL));
% $$$             ztmp = get_depths(hname,5,0,tind+ti-1);
% $$$             z_w(:,:,:,tind) = ztmp(xi:(xi+xL-1),yi:(yi+yL-1),zi:(zi+zL));
        end
    end
    eval(['var = ' GetVarStr ';']);
    %isopycnal plot:
    if (vtype == 4)
%Isopycnal:
        eval(['rho = (sw_dens0((' ...
                      'netcdf.getVar(' num2str(ncid) ',netcdf.inqVarID(' ...
                      num2str(ncid) ',''salt''),' Dom(ncid,'salt') ...
                  ',''double'')),(' ...
                      'netcdf.getVar(' num2str(ncid) ',netcdf.inqVarID(' ...
                      num2str(ncid) ',''temp''),' Dom(ncid,'temp') ...
                  ',''double'')))-1000);']);
% $$$         %Isothermal:
% $$$         eval(['rho = (netcdf.getVar(' num2str(ncid) ',netcdf.inqVarID(' ...
% $$$                       num2str(ncid) ',''temp''),' Dom(ncid,'temp') ...
% $$$                   ',''double''));']);

        NaNs = abs(rho(:,:,end,:))>100;
                    
        [X,Y] = meshgrid((1:xL),(1:yL));
        X = repmat(X',[1 1 zL]);
        Y = repmat(Y',[1 1 zL]);
        for t=1:tL
            FV = isosurface(X,Y,z_rho(:,:,:,t),rho(:,:,:,t),rhoISP, ...
                            var(:,:,:,t),'noshare');
            var(:,:,end,t) =  griddata(FV.vertices(:,1),FV.vertices(:,2), ...
                                     FV.facevertexcdata,X(:,:,1), ...
                                     Y(:,:,1));
        end
        var = var(:,:,end,:);
        var(NaNs) = NaN;
    end
    
    %unextend range if needed:
    if (uext == 1)
        var = var(2:(end-1),:,:,:);
    end
    if (vext == 1)
        var = var(:,2:(end-1),:,:);
    end
    if (zext == 1)
        var = var(:,:,2:(end-1),:);
    end
    if (textn == 1)
        if (vtype == 2)
            var(:,:,end+1) = NaN*var(:,:,end);
        else
            var(:,:,:,end+1) = NaN*var(:,:,:,end);
        end
    end
    netcdf.close(ncid);
    if (dname ~= 0)
        netcdf.close(ncidD);
    end
    if (exist(ename) == 2)
        netcdf.close(ncidE);
    end

    
function [OUT] = Drv(metric,Field,TYPE,varargin)
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
        metric = Int(metric,[1 1],[hgrdD 1]);
        
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
    OUT = Int(OUT,[hgrdD vgrdD],[1 1]);
end

function [OUT] = Int(Field,ingr,outgr)
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


function dstr = Dom(ncid,Vname)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function determines the domain string for the netcdf.getVar
% command given input ncid file identifier, variable name and vtype
% for output.
%
% INPUTS:
%
% ncid = id number of netcdf file.
%
% Vname = name of variable.
%
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 17/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------
    [hgrd,vgrd,Vdims] = grd_var(ncid,Vname);
if (Vdims == 2)
    switch hgrd
      case 1 % rho
        dstr = '[xi-1 yi-1],[xL yL]';
      case 2 % u
        dstr = '[xi-1 yi-1],[xL-1 yL]';
      case 3 % v
        dstr = '[xi-1 yi-1],[xL yL-1]';
      case 4 % psi
        dstr = '[xi-1 yi-1],[xL-1 yL-1]';
    end
elseif (Vdims == 3)
    switch hgrd
      case 1 % rho
        dstr = '[xi-1 yi-1 ti-1],[xL yL tL]';
      case 2 % u
        dstr = '[xi-1 yi-1 ti-1],[xL-1 yL tL]';
      case 3 % v
        dstr = '[xi-1 yi-1 ti-1],[xL yL-1 tL]';
      case 4 % psi
        dstr = '[xi-1 yi-1 ti-1],[xL-1 yL-1 tL]';
    end
elseif (Vdims == 4)
    switch hgrd
      case 1 % rho
        dstr = '[xi-1 yi-1 zi-1 ti-1],[xL yL ?? tL]';
      case 2 % u
        dstr = '[xi-1 yi-1 zi-1 ti-1],[xL-1 yL ?? tL]';
        dstr = '[xi yi-1 zi-1 ti-1],[xL-1 yL ?? tL]';

      case 3 % v
        dstr = '[xi-1 yi-1 zi-1 ti-1],[xL yL-1 ?? tL]';
        dstr = '[xi-1 yi zi-1 ti-1],[xL yL-1 ?? tL]';

      case 4 % psi
        dstr = '[xi-1 yi-1 zi-1 ti-1],[xL-1 yL-1 ?? tL]';
    end
    
    switch vgrd
      case 1 % rho
        dstr = strrep(dstr,'??','zL');
      case 2 % w
        dstr = strrep(dstr,'??','zL+1');
    end
end

end

function [hgrdnum,vgrdnum,NDimensions] = grd_var(ncid,varname)
%---------------------------------------------------------------------
%---------------------------------------------------------------------
%
% This function determines the grid numbers for the variable
% varname in the file fname.
%
% INPUTS:
%
% varname = name of variable string.
%
% fname = filenames
%
% OUTPUTS:
%
% hgrdnum = horizontal grid number 1(rho), 2(u), 3(v), 4(psi)
%
% vgrdnum = vertical grid number 1(rho) 2(w)
%
% NDimensions = number of dimensions
%---------------------------------------------------------------------
%
% Dependencies; none. 
%
% Ryans ROMS Matlab and netcdf Utilities 23/7/13
%
%---------------------------------------------------------------------
%---------------------------------------------------------------------
    if (strcmp(varname,'rho'))
        varname = 'temp';
    end

    vid = netcdf.inqVarID(ncid,varname);
    [tmp,tmp2,Dids,tmp3] = netcdf.inqVar(ncid,vid);
    NDimensions = length(Dids); % number of dimensions.
        
%Determine horizontal grid type:
    [dnh,tmp] = netcdf.inqDim(ncid,Dids(1));
    [dnh2,tmp] = netcdf.inqDim(ncid,Dids(2));

    if (dnh(end) == 'o')
        if (dnh2(end) == 'v')
            hgrdnum = 3;
        else
        hgrdnum = 1;
        end
    elseif (dnh(end) == 'u')
        hgrdnum = 2;
    elseif (dnh(end) == 'v')
        hgrdnum = 3;
    elseif (dnh(end) == 'i')
        hgrdnum = 4;
    end
    %Determine vertical grid type:
    vgrdnum = 1;
    if (NDimensions == 4)
        [dnv,tmp] = netcdf.inqDim(ncid,Dids(3));
        if (dnv(end) == 'w')
            vgrdnum = 2;
        end
    end
end



end