function Comm = init_geometry(GL,Nodes,indexord,thick)
% Comm = init_geometry(GL,Nodes,indexord,thick)
%
% INPUTS:
% GL(mu) = global lattice extension in direction mu;
% Nodes(mu) = number of computing nodes in which direction mu is distributed;
% indexord is any permutation of [1:D]. It defines completely which coordinates run slow or fast. 
%       The fastest running direction is indexord(1) and the slowest indexord(D).
% thick=depth of the border (allows the *direct* construction of observables of linear size 'thick').
%
% OUTPUTS:
% Comm is a structure with the following subfields.
% L(mu) = extension of the local sublattice in direction mu; i.e. L(mu)=GL(mu)/Nodes(mu);
% EL(mu) = extension in mu-direction of the  'extended lattice' = local lattice + borders.
% i.e. EL(mu) = L(mu)+2*thick, if direction mu is parallelized; 
%      EL(mu) = L(mu), if direction mu is not parallelized; 
% Ext_Ind(pnt) is the index (from 1 to EV=prod(EL)) in the extended lattice associated to the 
%       index pnt (from 1 to V=prod(V)) in the local sublattice. 
% Neigh(epnt,dir,ver) is the index in the extended lattice (from 1 to EV=prod(EL)) of the 
%       nearest Neighbours of epnt, when such neighbour also falls in the extended lattice 
%       (otherwise it is nan).
% Neigh_Node(mu,up/down) is the rank of the neighbouring node (to the one calling the function) 
%       in direction mu.
% BndInd is a structure that containes the indices of the borders distributed in a way which is 
%       useful for exchange_fields (see below for the precise organization).
%
%%%%%%% Note:
% dir is one of the D directions; ver=1(up),2(down);
%
% see also coord_basis, index2coord, coord2index, exchange_fields

L=GL./Nodes;
D=length(L);
dirs=[1:D];
par_dirs=find(Nodes>1);
EL=L; EL(par_dirs)=L(par_dirs)+2*thick; % Extended Lattice including borders where they are necessary.
V=prod(L);
EV=prod(EL);
cb=coord_basis(L,indexord);
ecb=coord_basis(EL,indexord);
viel_bein = eye(D,D);
sig=[1,-1];
if(any(thick>=L))
 exceed=find(any(thick>=L));
 sprintf('Borders exceeding the bulk in some directions: %d %d', exceed) 
end

%%% 1. Compute Ext_Ind(V)
for pnt=1:V
 coord = index2coord(pnt,cb,indexord);
 coord(par_dirs) = coord(par_dirs)+thick; % (x_mu) -> (x_mu + thick) if mu is parallelized.
 Ext_Ind(pnt) = coord2index(coord,ecb);
end

%%% 2. Compute Neigh(EV,D,2)
for epnt=1:EV
  ecoord=index2coord(epnt,ecb,indexord);
  for mu=1:D
   for ver=1:2
    ncoord=ecoord+sig(ver)*viel_bein(mu,:);
    if(all(ncoord(par_dirs)<EL(par_dirs)) && all(ncoord(par_dirs)>=0))% if the nearest-neighbouring point of epnt
								      % (which may already be in the border)
								      % is still in the extended lattice.
     Neigh(epnt,mu,ver)=coord2index(mod(ncoord,EL),ecb); % define Neigh
    else
     Neigh(epnt,mu,ver)=nan;
    end
   end
  end
end

%%%%%%%%%%%%%%%%% The following is needed only in the parallel case.  %%%%%%%%%%%%%%%%%%%%%%%%%%
if ( any(L ~= EL) )
%%%% 3. Compute Neigh_Node(D,up/dn), which returns the rank of the nodes neighbouring to this one.

 NPD=length(par_dirs);
 send=1; recv=2; up=1; dn=2;
 Ncb=coord_basis(Nodes,indexord);

 % MPI_Init;
 % [info mpi_rank] =  MPI_Comm_rank(MPI_COMM_WORLD);
 mpi_rank =0;
 rank=mpi_rank+1; % check that rank returned by MPI_Comm_rank actually goes from 0 to num_proc-1 ???

 for mu=1:D
  if (EL(mu)==L(mu))
    Neigh_Node(mu,up)=rank;
    Neigh_Node(mu,dn)=rank;
  else
    Neigh_Node(mu,up)=coord2index(mod(index2coord(rank,Ncb,indexord)+viel_bein(mu,:),Nodes),Ncb);
    Neigh_Node(mu,dn)=coord2index(mod(index2coord(rank,Ncb,indexord)-viel_bein(mu,:),Nodes),Ncb);
  end
 end

%%%% 4. Define the struct BndInd  containing the indices of the fields to be communicated

 for mu=par_dirs; % loop on the parallelized dimensions
  for en=1:NPD
   eval(['list_dn_recv_' num2str(en) '= [];' ]); % set temporary variables list_*_*_* to []
   eval(['list_up_recv_' num2str(en) '= [];' ]); 
   eval(['list_dn_send_' num2str(en) '= [];' ]);
   eval(['list_up_send_' num2str(en) '= [];' ]);
  end
  other_par_dirs=setdiff(par_dirs,mu);
  oned=ones(1,length(other_par_dirs));
  for epnt=1:EV
   ecoord=index2coord(epnt,ecb,indexord);
   other_par_coord=ecoord(other_par_dirs); 
   edg_num=...                             
       length(find((ecoord(other_par_dirs) < thick*oned) | ...
		   (ecoord(other_par_dirs) >= EL(other_par_dirs)-thick*oned)))+1;
   % edg_num is the number of components in ecoord that live in the  border of the extended lattice. 
   % e.g. corners of a cube have edg_num=3. Note that I need at least edg_num communications steps to 
   % update them; that's why I divide BndInd according to edg_num.
   if (ecoord(mu) < thick)
    eval(['list_dn_recv_' num2str(edg_num) '= [list_dn_recv_' num2str(edg_num) '; epnt];' ])
   elseif ((thick <= ecoord(mu)) && (ecoord(mu) < 2*thick))
    eval(['list_dn_send_' num2str(edg_num) '= [list_dn_send_' num2str(edg_num) '; epnt];' ])
   elseif ((EL(mu)-2*thick <= ecoord(mu)) && (ecoord(mu) < EL(mu)-thick))
    eval(['list_up_send_' num2str(edg_num) '= [list_up_send_' num2str(edg_num) '; epnt];' ])
   elseif (ecoord(mu) >= EL(mu)-thick)
    eval(['list_up_recv_' num2str(edg_num) '= [list_up_recv_' num2str(edg_num) '; epnt];' ])
   end
  end % epnt= 1:EV

  for en=1:NPD
   eval(['BndInd.en_' num2str(en) '_mu_' num2str(mu) '(:,dn,recv)= list_dn_recv_' num2str(en) ';'])
   eval(['BndInd.en_' num2str(en) '_mu_' num2str(mu) '(:,dn,send)= list_dn_send_' num2str(en) ';'])
   eval(['BndInd.en_' num2str(en) '_mu_' num2str(mu) '(:,up,send)= list_up_send_' num2str(en) ';'])
   eval(['BndInd.en_' num2str(en) '_mu_' num2str(mu) '(:,up,recv)= list_up_recv_' num2str(en) ';'])
  end
 end % imu=1:NPD

% Note. BndInd is a struct, because matlab does not allow dimension mismatch in a single array (but the number 
% of indices associated to different directions are necessarely different). Moreover, struct field names cannot  
% be numerical variables. Hence, I had to use the ugly eval()'s above.

end % (EL~=L)

% 5. define Comm structure in order to collect all variables related to neighbour addressing and communications

Comm.l=L;
Comm.el=EL;
Comm.pd=par_dirs;
Comm.neigh=Neigh;
Comm.extind=Ext_Ind;
Comm.flag=0;
if (EL~=L)
 Comm.neighnode=Neigh_Node;
 Comm.bndind=BndInd;
 Comm.mpirank=mpi_rank; % do I need also MPI_COMM_WORLD ?
 Comm.flag=1;
end
