function [psi] = init_spinor_field(Nc, Comm, init_flag, norm,filename)
% init_flag      -1 (readin from filaname) 
%                 0 (filled with zeros)
%                 1 (gaussian and globally normalized to norm)

L=Comm.l;
EL=Comm.el;
Ext_Ind=Comm.extind;
V=prod(L);
EV=prod(EL)
D=length(L);
S=2^(floor(D/2));
psi=zeros(Nc,S,EV);
if (init_flag == 1)
  psi(:,:,Ext_Ind) = (randn(Nc,S,V) + i * randn(Nc,S,V));
  psi(:,:,Ext_Ind) = psi(:,:,Ext_Ind)/sqrt(V*S*Nc*2);
elseif(init_flag==-1)
  psi=load(filename);
end
