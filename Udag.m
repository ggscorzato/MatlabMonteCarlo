function Ud = Udag(U,Comm)
% Ud = Udag(U,L,Neigh,Ext_Ind)
% Compute U^dag in order to optimize a version of the fermion matrix.

UP=1;
L=Comm.l;
D=length(L); V=prod(L);
for mu=1:D
 n_up = Comm.neigh(Comm.extind,mu,UP);
 Ud(:,:,mu,Comm.extind) = U(:,:,mu,n_up);
end
Ud = conj(permute(Ud,[2 1 3 4]));
