function psi1 = fermion_matrix(psi0,U,Ud,Comm,Nc,mass,gi,gv)
% usage:  psi1 = fermion_matrix(psi0,U,Ud,Comm,Nc,mass,gi,gv)
% where:  Comm=communicator (see init_geometry);
%         Nc = number of colors (SU(Nc) gauge group in fund. rep);
%         U= gauge configuration. It must be a structure U(D,V).m(Nc,Nc);
%         Ud = U^dag obtained with Udag.m (borders are not needed).
%         psi0= spinorial vector psi0(Nc,S,V);
%         S= dim of Dirac spin rep. (=2^(floor(D/2));

rW=1;
EL=Comm.el;
Ext_Ind=Comm.extind;
EV=prod(EL);
D=length(EL);
S=2^(floor(D/2));
UP=1;
DN=2;


psi1(:,:,Ext_Ind)=(D*rW + mass) * psi0(:,:,Ext_Ind);

for mu=1:D
  for s=1:S
    p_up=Neigh(Ext_Ind,mu,UP);
    p_dn=Neigh(Ext_Ind,mu,DN);
    Xid(:,s,Ext_Ind,mu) = rW * psi0(:,s,p_dn) + gv(mu,s) * psi0(:,gi(mu,s),p_dn);
    Xiu(:,s,Ext_Ind,mu) = rW * psi0(:,s,p_up) - gv(mu,s) * psi0(:,gi(mu,s),p_up);
  end
end

for s=1:S
  for c=1:Nc
   tens(:,Ext_Ind,:,c,s) =  (squeeze(Ud(c,:,:,Ext_Ind)) .* squeeze(Xid(:,s,Ext_Ind,:)) + ...
			     squeeze(U(c,:,:,Ext_Ind)) .* squeeze(Xiu(:,s,Ext_Ind,:)));
  end
end

for s=1:S
  for c=1:Nc
    tens1(Ext_Ind,c,s) = squeeze(sum(sum(tens(:,Ext_Ind,:,c,s),3),1));
  end
end
psi1(:,:,Ext_Ind) = psi1(:,:,Ext_Ind) + permute(tens1(Ext_Ind,:,:),[2,3,1]);

% communications (no need of edges)
if (Comm.flag ==1)
  psi1=reshape(exchange_fields(reshape(psi1,Nc*S,EV),Comm,1),Nc,S,EV);
end