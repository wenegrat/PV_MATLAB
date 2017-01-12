function OUT = DrvS(metric,zmetric, Field,TYPE)
if (TYPE == 'z') 
    OUT = DrvROMS(metric, Field, TYPE);
else
    dFdx = DrvROMS(metric, Field, TYPE);
    dZdx = DrvROMS(metric, zmetric, TYPE);
    dFdz = DrvROMS(zmetric, Field, 'z');
    OUT = dFdx - dFdz.*dZdx;
end

end