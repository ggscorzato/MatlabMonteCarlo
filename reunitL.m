function u = reunitL(v)
% u = reunitL(v)
% reunitarizaion according to Liang et al. Phys.Lett.B 307 (1993) 375
Nc=size(v,1);
v1=v';
v1v=inv(sqrtm(v1*v));
d=det(inv(v) * v1)^(1/6);
u=v * v1v * d;
