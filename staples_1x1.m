function stap = staples_1x1(U,Comm,needforlarger)
% stap = staples_1x1(U,Comm,needforlarger)
%
% compute the staples 1x1 for the configuration U and communication structure Comm.
% If needforlarger =1, the result is saved in a vector with indices:
% stap(Nc,Nc,direction of the hingeing link, direction of the staple, up/down, extended volume)
% If needforlarger =0, the staples are averaged over the 2*D possibilities and results are saved in:
% stap(Nc,Nc,direction of the hingeing link,internal volume)


up=1;dn=2;
dirs=[1:length(Comm.l)];
EV=prod(Comm.el);
V=prod(Comm.l);
Nc=size(U,1);
stap=zeros(Nc,Nc,D,D,2,EV);

for mu= dirs
 for nu = setdiff(dirs,mu)
  for pnt =1:V
   epnt=Comm.extind(pnt);
   stap(:,:,mu,nu,up,epnt)=...
       U(:,:,nu,epnt)*...
       U(:,:,mu,Comm.neigh(epnt,nu,up))*...
       U(:,:,nu,Comm.neigh(epnt,mu,up))';
   stap(:,:,mu,nu,dn,epnt)=...
       U(:,:,nu,Comm.neigh(epnt,nu,dn))' *...
       U(:,:,mu,Comm.neigh(epnt,nu,dn))*...
       U(:,:,nu,Comm.neigh(Comm.neigh(epnt,nu,dn),mu,up));

  end
 end
end

if ( needforlarger == 1 )
 stap = reshape(stap,Nc*Nc*D,D,2,EV);
 stap = exchange_staples(stap,Comm);
 stap = reshape(stap,Nc,Nc,D,D,2,EV);
elseif ( needforlarger == 0 )
 stap=squeeze(mean(mean(stap,4),5)); % mean over the 2*D directions
end
