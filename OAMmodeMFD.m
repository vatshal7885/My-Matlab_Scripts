function [ex ey] = OAMmodeMFD(TH,R,mfd,l,m)
mfd = mfd/2;
mfd = mfd/sqrt(2);
V = 1/(mfd.^2);

LG = LaguerrePoly3((m-1),l);

F_lm = R.^l.*polyval(LG,V.*R.^2).*exp(-0.5*V.*R.^2);
l_TH = l.*TH;
Ex = F_lm.*cos(l_TH);
Ey = F_lm.*sin(l_TH);


ex = Ex+1i.*Ey;
ey = Ex-1i.*Ey;

ex = ex./sqrt(sum(sum(abs(ex).^2)));
if (l>0)
    ey = ey./sqrt(sum(sum(abs(ey).^2)));
end
% Calculates the Laguerre polynomials
function [coeff] = LaguerrePoly3(n,alpha)
coeff = zeros(1,n);
FACTORIAL_1 = factorial((n+alpha));
for m=0:n
    bino = FACTORIAL_1/(factorial((n-m))*factorial((n+alpha)-(n-m)));
    coeff((n+1)-m)=bino*(-1).^m/(factorial(m));
end