function U = init_gauge_fields(Nc, init_flag,alpha,filename,Comm)
% U = init_gauge_fields(Nc, init_flag,alpha,filename,Comm)
%
% Initialize a gauge configuration with gauge group SU(Nc)
% init_flag     1 (random with spread alpha, cold for alpha=0)
% init_flag    -1 (readin filename)
% Comm is the communicator structure (see init_geometry.m)
% Indices are ordered according to: U(Nc,Nc,EV,D).

Neigh=Comm.neigh;
Ext_Ind=Comm.extind;
L=Comm.l;
EL=Comm.el;
V=prod(L);
EV=prod(EL);
D=length(L);
ge=sun_gen(Nc,1);

if (init_flag == 1)
  % generate random gauge field
  for pnt=1:V
    for mu=1:D
      a=0;
      for c=1:Nc^2-1
	a = a + ge(:,:,c) * ((rand-0.5)*alpha); 
      end
      U(:,:,mu,Ext_Ind(pnt)) = expm(a); 
    end
  end
elseif(init_flag ==-1)
  % read gauge field from disk.
  conf=load(filename);
  U=conf.U;
  clear conf
end

% update boundaries
max_num_edges=2; % unless you need some stange higher dim observable
U = reshape(exchange_fields(reshape(U,Nc*Nc*D,EV),Comm,max_num_edges),Nc,Nc,D,EV);