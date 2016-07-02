function ddx = DPeriodic(input, dx, dy, dir)
if strcmp(dir, 'x')
    ddx = DxPeriodic(input, dx);
end
if strcmp(dir, 'y')
    ddx = DyPeriodic(input, dy);
end
end

function didx = DxPeriodic(input, dx)
%Assumes inputs are on rho points.
didx = Drv(dx, input, 'x');

endcap = (input(1,:,:,:) - input(end,:,:,:))./dx; % sits on U point = 1.
rcap = (input(2,:,:,:) - input(1,:,:,:))./dx; % sits on U point = 2;
lcap = (input(end,:,:,:) - input(end-1,:,:,:))./dx; %sits on U point end-1

didx(1,:,:,:) = 0.5.*(rcap + endcap); %average back to rho=1 point.
didx(end,:,:,:) = 0.5*(endcap+lcap); %average back to rho=end point
end