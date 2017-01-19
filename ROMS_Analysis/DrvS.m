function OUT = DrvS(metric,zmetric, Field,TYPE, varargin)
if (TYPE == 'z') 
    OUT = DrvROMS(metric, Field, TYPE, varargin{:});
else
    dFdx = DrvROMS(metric, Field, TYPE, varargin{:});
    dZdx = DrvROMS(metric, zmetric, TYPE, 1);
    dFdz = DrvROMS(zmetric, Field, 'z', varargin{:});
    OUT = dFdx - dFdz.*dZdx;
end

end