function U = hyp_smearing(U,Nc,alpha,Comm)
%  U = hyp_smearing(U,Nc,apha,Comm)
%
%  see Hasenfratz Knechtli Phys.Rev.D64 034504
%  generalized to D dimensions
%

L=Comm.l;
V=prod(L);
EL=Comm.el;
EV=prod(EL);
D=length(L);
dirs=[1:D];
if (isempty(alpha))
 alpha=zeros(D-1);
 alpha=[0.75 0.6 0.3];
end

jmax=0;
for nex = D-2:-1:0
    jmax=jmax+D*nchoosek(D-1,nex);
end  
extsizeu=size(U);
extsizeu(3)=jmax;
UU= zeros(extsizeu);
UU(:,:,1:D,:)=U;

for nex = D-2:-1:0  % loop on the smearing iterations: U->\tildaV->\hatV->V. 
		    % nex= number of dirs excluded from smearing besides mu.
  thesej=[];
  for mu = 1:D        % loop on direction of the new link mu
    orth = dirs(find(dirs)~=mu); % = list of all non-mu directions
    nurho=nchoosek(orth,nex); % = list of possible choices of nex non-mu directions
    nt=nchoosek(D-1,nex); % = number of possible choices of nex non-mu directions
    for t = 1:nt           % loop on the possible choice of nex non-mu directions
      j=V_index(D,mu,nurho(t,:)); 
      eta=setdiff(orth,nurho(t,:));
      thesej=[thesej,j];
      for h = eta % loop over the remaining (D-1-nex) directions (over which the sums are performed).
	j1(h)=V_index(D,mu,sort([h  nurho(t,:)])); 
	j2(h)=V_index(D,h ,sort([mu nurho(t,:)]));
      end
      for epnt = Comm.extind      % loop on space-time
	stap=zeros(Nc,Nc);
	for h = eta      % sum on the remaining (D-1-nex) (positive and negative) directions 
	  stap=stap+...
	       UU(:,:,j2(h),epnt)*...
	       UU(:,:,j1(h),Comm.neigh(epnt,h,1))*...
	       UU(:,:,j2(h),Comm.neigh(epnt,mu,1))'+...
	       UU(:,:,j2(h),Comm.neigh(epnt,h,2))' *...
	       UU(:,:,j1(h),Comm.neigh(epnt,h,2))*...
	       UU(:,:,j2(h),Comm.neigh(Comm.neigh(epnt,h,2),mu,1));
	end % h
	UU(:,:,j,epnt)=reunit((1-alpha(nex+1))*U(:,:,mu,epnt) + ...
				     (alpha(nex+1)/(2*length(eta)))*stap);
      end % epnt
    end % t
  end % mu
  ltj=length(thesej);
  UU(:,:,thesej,:)=reshape(exchange_fields(reshape(UU(:,:,thesej,:),Nc*Nc*ltj,EV),Comm,ltj),Nc,Nc,ltj,EV);
end % nex
U=UU(:,:,1:D,:);

function j = V_index(D,mu,v)
% function that maps the indices {mu;nu,rho,..} into a sequential index j
% (keeping all indices would waste a lot of storage since only mu~=nu~=rho~=... are needed).
sv=sort(v);
lv=length(sv);
if (lv==D-1) % U_mu (the starting links)
  j=mu;
elseif(lv==0) % V_mu (the final links, which overwrite the starting links)
  j=mu;
else % the indices of the others are more complicated
 j=0;
 for ll=lv+1:D-1
   j=j+D*nchoosek(D-1,ll);
 end
 dirs=[1:D];
 nt=nchoosek(D-1,lv);
 orth = dirs(find(dirs)~=mu);
 p=nchoosek(orth,lv);
 j=j+(mu-1)*nchoosek(D-1,lv);
 j=j+find(prod(double(p==ones(nt,1)*sv),2));
end
% for D=4 we have:
% mu v         j
% 1 []/[2 3 4] 1   U_mu and V_mu
% 2 []/[1 3 4] 2
% 3 []/[1 2 4] 3
% 4 []/[1 2 3] 4
% 1 [2 3]      5   Vbar_mu;[v]
% 1 [2 4]      6
% 1 [3 4]      7
% 2 [1 3]      8
% 2 [1 4]      9
% 2 [3 4]      10
% 3 [1 2]      11
% 3 [1 4]      12
% 3 [2 4]      13
% 4 [1 2]      14
% 4 [1 3]      15
% 4 [2 3]      16
% 1 [2]        17  Vtild_mu;[v]
% 1 [3]        18
% 1 [4]        19
% 2 [1]        20
% 2 [3]        21
% 2 [4]        22
% 3 [1]        23
% 3 [2]        24
% 3 [4]        25
% 4 [1]        26
% 4 [2]        27
% 4 [3]        28
