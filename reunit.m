function u = reunit(v)
% u = reunit(v)
% unitaize a square matrix by anti-symmetrizing the logarithm
Nc=size(v,1);
a=logm(v);
a=0.5*(a-a');
a=a-eye(Nc)*trace(a)/Nc;
u=expm(a);
