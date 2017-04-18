function ddx = DPeriodic(input, dl, dir)
if strcmp(dir, 'x')
    ddx = DxPeriodic(input, dl);
end
if strcmp(dir, 'y')
    ddx = DyPeriodic(input, dl);
end
end

function didx = DxPeriodic(input, dx)
%Assumes inputs are on rho points.
nx= size(input, 1);
if nx>1
didx = Drv(dx, input, 'x');

endcap = (input(1,:,:,:) - input(end,:,:,:))./dx; % sits on U point = 1.
rcap = (input(2,:,:,:) - input(1,:,:,:))./dx; % sits on U point = 2;
lcap = (input(end,:,:,:) - input(end-1,:,:,:))./dx; %sits on U point end-1

didx(1,:,:,:) = 0.5.*(rcap + endcap); %average back to rho=1 point.
didx(end,:,:,:) = 0.5*(endcap+lcap); %average back to rho=end point

else
    didx = 0.*input;
end
end

function didx = DyPeriodic(input, dy)
didx = Drv(dy, input, 'y');

endcap = (input(:,1,:,:) - input(:,end,:,:))./dy;%sits at V point = 1.
scap = (input(:,2,:,:) - input(:,1,:,:))./dy; % V point = 2
ncap = (input(:,end,:,:) - input(:,end-1,:,:))./dy; % vpoint end-1
didx(:,1,:,:) = 0.5*(scap +endcap);
didx(:,end,:,:) = 0.5*(ncap+endcap);
end