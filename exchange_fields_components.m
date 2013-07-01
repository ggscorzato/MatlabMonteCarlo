function FF = exchange_fields_components(FF,Comm,IndexComp)
% FF = exchange_fields_components(FF,Comm,IndexComp)
%
% takes a field FF defined on each node and copies the boundaries of each node onto the mirror boundaries of
% neighbours nodes. 
%
% This function copies only the components defined below.
%
% FF is assumed to be an 2 dimensional array. The first dimension can be of any size and it is copied as it is 
% (here we think of spin, color, lorentz, replica, etc. A reshape/permute is tipically necessary). 
% The second dimension represents the space-time index.
%
% Comm is the communicator structure defined in init_geometry.
% IndexComp selacts the components that should be sent
%
% see also init_geometry, exchange_fields, index_staples
% 

send=1; recv=2; up=1; dn=2;
size_of_ff=64 * size(FF,1); %% ??? check with mpitb if size_of_double (=64) needs to be specified.

en=1;
for mu=Comm.pd
  eval(['bndi=Comm.BndInd.en_' num2str(en) '_mu_' num2str(mu) ';'])
  siz=length(bndi(:,1,1));
  %    send mu+ 
  MPI_Send(FF(:,IndexComp(mu,up,send),bndi(:,up,send)),siz,size_of_ff,Comm.NeighNode(this_node,mu,up),...
	   98,MPI_COMM_WORLD);
  %    recv mu-
  MPI_Recv(FF(:,IndexComp(mu,dn,recv),bndi(:,dn,recv)),siz,size_of_ff,Comm.NeighNode(this_node,mu,dn),...
	   98,MPI_COMM_WORLD);
  %    send mu-
  MPI_Send(FF(:,IndexComp(mu,dn,send),bndi(:,dn,send)),siz,size_of_ff,Comm.NeighNode(this_node,mu,dn),...
	   99,MPI_COMM_WORLD);
  %    recv mu+
  MPI_Recv(FF(:,IndexComp(mu,up,recv),bndi(:,up,recv)),siz,size_of_ff,Comm.NeighNode(this_node,mu,up),...
	   99,MPI_COMM_WORLD);
end
