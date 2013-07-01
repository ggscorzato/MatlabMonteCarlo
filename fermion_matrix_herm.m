function psi1 = fermion_matrix_herm(psi0,U,Ud,Comm,Nc,mass,gi,gv)
% usage:  psi1 = fermion_matrix_herm(psi0,U,Ud,Comm,Nc,mass,gi,gv);
% where:  L = [LX,LY,LZ,T...];
%         LX = number of point in direction X, LY= ...
%         D=length(L);
%         Nc = number of colors (SU(Nc) gauge group in fund. rep);
%         U= gauge configuration.
%         psi0= spinorial vector psi0(Nc,S,V);
%         S= dim of Dirac spin rep. (=2^(floor(D/2));

L=Comm.l;
D=length(L);
S=2^(floor(D/2));

psi_tmp=fermion_matrix(psi0,U,Ud,Comm,Nc,mass,gi,gv);
for s=1:S
 psi1(:,s,:,:)=gv(D+1,s) * psi_tmp(:,gi(D+1,s),:,:); 
end