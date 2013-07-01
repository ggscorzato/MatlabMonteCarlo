function U_BS = border_staples(U,Comm)
% U_BS = border_staples(U,Comm)
%
% Compute staples hanging on border links but extending on the interior.
% These are needed to avoid using a border of depth=2 when computing 2x1 wilson loops.

send=1; recv=2; up=1; dn=2; 
L=Comm.l; D=length(L); all_dirs=[1:D];
for mu = Comm.pd
 for bord= [dn up]
  eval(['ind=Comm.bndind.en_1_mu_' num2str(mu) '(:,bord,send);'])
  for nu=setdiff(all_dirs,mu)
    %%% SEGNO
  end
 end
end
