function KE = calculateKE(U, V, W,X, Y, Z,T)

x = linspace(X(1), X(end), 100);
y = linspace(Y(1), Y(end), 100);
z = linspace(Z(1), Z(end), 100);

Ured = interpn(X, Y, Z, T, U, x, y, z, T);
Vred = interpn(X, Y, Z, T, V, x, y, z, T);
Wred = interpn(X, Y, Z, T, W, x, y, z, T);

KE = 1035.*(Ured.^2 +Vred.^2 + Wred.^2)./2;

KE = squeeze(nansum(nansum(trapz(z, KE, 3))));

end