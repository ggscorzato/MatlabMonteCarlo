function FF = exchange_fields(FF,Comm,MaxNumEdg)
% FF = exchange_fields(FF,Comm,MaxNumEdg)
%
% takes a field FF defined on each node and copies the boundaries of each node onto the mirror boundaries of
% neighbours nodes. 
%
% FF is assumed to be an 2 dimensional array. The first dimension can be of any size and it is copied as it is 
% (here we think of spin, color, lorentz, replica, etc. A reshape/permute is tipically necessary). 
% The second dimension represents the space-time index.
%
% Comm is the communicator structure defined in init_geometry.
%
% MaxNumEdg is the higher level of corner that we need to communicate. E.g. For Wilson fermions it is 1; 
% for Wilson loops it is 2; for D-dim hyp-smearing it is D.
% 
% see also init_geometry
% 
% ??? optimization
% matlab advertizes its  "copy on write" that should avoid the duplication of unmodified variables (FF in this 
% case). but I am skeptic in this case because: 1. I modify part of an array and 2. I copy back to FF.
% I guess I will have to define global fields and distinct exchange_** functions ... (give me pointers !!!!)

NPD=length(Comm.pd);
send=1; recv=2; up=1; dn=2;
size_of_ff=64 * size(FF,1); %% ??? check with mpitb if size_of_double (=64) needs to be specified.

for en=1:min(NPD,MaxNumEdg)
 for mu=Comm.pd
  eval(['bndi=Comm.BndInd.en_' num2str(en) '_mu_' num2str(mu) ';'])
  siz=length(bndi(:,1,1));
  %    send mu+ 
  MPI_Send(FF(:,bndi(:,up,send)),siz,size_of_ff,Comm.NeighNode(this_node,mu,up),98,MPI_COMM_WORLD);
  %    recv mu-
  MPI_Recv(FF(:,bndi(:,dn,recv)),siz,size_of_ff,Comm.NeighNode(this_node,mu,dn),98,MPI_COMM_WORLD);
  
  %    send mu-
  MPI_Send(FF(:,bndi(:,dn,send)),siz,size_of_ff,Comm.NeighNode(this_node,mu,dn),99,MPI_COMM_WORLD);
  %    recv mu+
  MPI_Recv(FF(:,bndi(:,up,recv)),siz,size_of_ff,Comm.NeighNode(this_node,mu,up),99,MPI_COMM_WORLD);
  end
end
