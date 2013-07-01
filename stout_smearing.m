function U = stout_smearing(U,Nc,Comm)
% U = stout_smearing(U,Nc,Comm)
% 
%  see Morningstar and Peardon hep-lat/0311018
%  generalized to D dimensions

up=1;dn=2;
L=Comm.l;
V=prod(L);
D=length(L);
rho=0.1;
rhom=rho * ones(4,4);

alld=[1:D];

for epnt = Comm.extind      % loop on space-time
  for mu = 1:D              % loop on direction of the new link mu
    C_mu=zeros(Nc,Nc);
    orth = alld(find(alld)~=mu); % = list of all non-mu directions
    for nu = orth       % loop on all the  non-mu directions
      C_mu = C_mu + ...
	     rhom(mu,nu) * (...
		 U(:,:,nu,epnt)*...
		 U(:,:,mu,Comm.neigh(epnt,nu,up))*...
		 U(:,:,nu,Comm.neigh(epnt,mu,up))'+...
		 U(:,:,nu,Comm.neigh(epnt,nu,dn))' *...
		 U(:,:,mu,Comm.neigh(epnt,nu,dn))*...
		 U(:,:,nu,Comm.neigh(Comm.neigh(epnt,nu,dn),mu,up)) );
    end
    Omega = C_mu * U(:,:,mu,epnt)';
    dOmega= Omega' - Omega;
    Qmu = 0.5 * dOmega - 0.5/Nc * trace(dOmega) *eye(Nc);
    UU(:,:,mu,epnt)=expm(Qmu) * U(:,:,mu,epnt);

  end
end

U=reshape(exchange_fields(reshape(UU,Nc*Nc*D,EV),Comm,1),Nc,Nc,D,EV)
