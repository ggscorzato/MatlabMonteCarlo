function wloop = wilson_loop(U,l1,dirs1,l2,dirs2,Comm)
% wloop = wilson_loop(U,l1,dirs1,l2,dirs2,Comm)
%
% Compute Wilson loop (= 1/V sum_x 1/n_dirs sum_{mu\in dirs1,nu \in dirs2} Re Tr[U_path(x,l1,mu,l2,nu)] ) 
% where path(x,l1,mu,l2,nu) is a rectangular path of lengths l1 in direction mu and l2 in direction nu, 
% starting from (its lowest) corner x.
% dirs1 and dirs2 must be row vectors with a list of directions. The directions contained in dirs1 and dirs2 
% must either the same set or be disjoint sets.
%
% see also polyakov_loop

%%%% 1. Define ll(dirs) and len_per_dir in order to deal with different cases uniformly.
V=prod(Comm.l);
EV=prod(Comm.el);
i1=1;i2=1;
if(dirs1==dirs2) % CASE A: square loops, involving the same set of directions (e.g. space-like or any squares).
 len_per_dir=1; % In this case I need to store only one Wilson Line per direction
 dirs=dirs1;
 ll(dirs)=max(l1,l2);
 if (l1<l2) % CASE B: rectangles build with a single set of directions (e.g. space-like or any rectangles).
  len_per_dir=2; % In this case I need to store 2 Wilson Lines per direction
  i1=2; i2=1;  % the first is always the shorter
 elseif(l2<l1)
  len_per_dir=2; 
  i1=1; i2=2;  % the first is always the shorter
 end
elseif(isempty(intersect(dirs1,dirs2))) % CASE C: rectangles or squares out of non-overlapping sets of 
					% directions (e.g. time-like loops).
 len_per_dir=1; % In this case I need to store only one Wilson Line per direction
 dirs=union(dirs1,dirs2);
 ll(dirs1)=l1;
 ll(dirs2)=l2;
else  % CASE D: all other cases can be build out of the previous A, B, C.
 sprintf('The row vectors dirs1 and dirs2 must be either the same set or disjoint sets. See help wilson_loop')
end

%%%%% 2. Construct all the needed Wilson Lines and leave them on the leftmost (lowest) end: 'outward stars'.
tmp=zeros(Nc,Nc,V);
Wline=zeros(Nc,Nc,EV,length(dirs),len_per_dir);

for imu=1:length(dirs)
 mu=dirs(imu);
 P_dn=Comm.neigh(Comm.extind,mu,2);
 Wline(:,:,:,imu,1)=squeeze(U(:,:,mu,:)); % First link of Wilson line in the mu direction
 if( (len_per_dir == 2) && (1==min(l1,l2)) ) % If you need more lengths for the same direction (case B),
					     % and if the shortest has length=1.
  Wline(:,:,:,imu,2)=Wline(:,:,:,imu,1);
  Wline(:,:,:,imu,2)=reshape(exchange_fields(reshape(Wline(:,:,:,imu,2),Nc*Nc,EV),Comm,1),Nc,Nc,EV); 
 end
 for m=2:ll(mu)
   for pnt=1:V
     % The Wilson line grows in the backward direction 
     % (note that tmp does not need to be defined on the boundaries).     
     tmp(:,:,pnt) = U(:,:,mu,P_dn(pnt)) * Wline(:,:,Comm.extind(pnt),imu,1); 
   end
   Wline(:,:,P_dn,imu,1)=tmp(:,:,:); % Wline lives at the lowest end, where the Wilson line grows.
				     % update and communicate borders of Wline
   Wline(:,:,:,imu,1)=reshape(exchange_fields(reshape(Wline(:,:,:,imu,1),Nc*Nc,EV),Comm,1),Nc,Nc,EV); 

   if( (len_per_dir == 2) && (m==min(l1,l2)) ) % If you need more lengths for the same direction (case B),
					      % and if you are done with the shortest.
    Wline(:,:,:,imu,2)=Wline(:,:,:,imu,1);
    Wline(:,:,:,imu,2)=reshape(exchange_fields(reshape(Wline(:,:,:,imu,2),Nc*Nc,EV),Comm,1),Nc,Nc,EV); 
   end
 end % m=1:ll(mu)
end % imu=1:length(dirs)
clear tmp

%%%% 3. Copy the Wilson Lines where they are needed and take the traces.

etmp=zeros(Nc,Nc,EV);
for imu1=1:length(dirs1)
 mu1=dirs1(imu1);
 P_up1=Comm.neigh(Comm.extind,mu1,1);
 for imu2=1:length(dirs2)
  if((dirs1~=dirs2) || (imu2>imu1) || (l1~=l2))
   mu2=dirs2(imu2);
   P_up2=Comm.neigh(Comm.extind,mu2,1);

   % bring the mu2 wilson line back of l1 steps
   etmp(:,:,Comm.extind)=Wline(:,:,Comm.extind,imu2,i2);
   for m=1:l1
    etmp(:,:,Comm.extind)=etmp(:,:,P_up1);
    etmp=reshape(exchange_fields(reshape(etmp,Nc*Nc,EV),Comm,1),Nc,Nc,EV);
   end

   % form the staple mu2 - mu1 - mu2
   for pnt=1:V
    epnt=Comm.extind(pnt);
    WL(:,:,pnt) = etmp(:,:,epnt) * Wline(:,:,epnt,imu1,i1) * Wline(:,:,epnt,imu2,i2)';
   end

   % bring the mu1 wilson line back of l2 steps
   etmp(:,:,Comm.extind)=Wline(:,:,Comm.extind,imu1,i1);
   for m=1:l2
    etmp(:,:,Comm.extind)=extind(:,:,P_up2);
    etmp=reshape(exchange_fields(reshape(etmp,Nc*Nc,EV),Comm,1),Nc,Nc,EV);
   end

   % close the staple with the last mu1
   for pnt=1:V
    epnt=Comm.extind(pnt);
    WL(:,:,pnt) = etmp(:,:,epnt)' * WL(:,:,pnt);
   end

   % Take the traces
   wloop(imu1,imu2)=mean(trace(WL(:,:,:)),3);

  end % if ( () || () || ())
 end % imu2
end % imu1
