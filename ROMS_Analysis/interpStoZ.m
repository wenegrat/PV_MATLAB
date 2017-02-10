function fieldI = interpStoZ(field, zm, z)
field = permute(field, [3 1 2]);
z = permute(z, [3 1 2]);

fieldI = interp1(z, field, zm);
fieldI = permute(fieldI, [2 3 1]);
end