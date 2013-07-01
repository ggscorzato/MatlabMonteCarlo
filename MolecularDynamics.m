function U = MolecularDynamics(U,@Mh,beta,ge,Nc,dt,N_MD,Comm)
% U = MoleculerDynamics(U,@Mh,beta,ge,Nc,dt,N_MD,Comm)
%

L=Comm.l;
D=length(L);
V=prod(L);
S=2^(floor(D/2));
EV=prod(Comm.el);
tol=1e-6;
maxit=1e4;
comloc=Comm;
Nadj=Nc^2-1;
epsilon=1e-6;
UP=1;
Ud=Udag(U,Comm);

for a=1:Nadj
  uge(:,:,a)=expm(i * epsilon * ge(:,:,a));
end

% 1. Pseudofermions
R_field = (randn(Nc,S,EV) + i * randn(Nc,S,EV))/sqrt(2); % check normalization ?!?
Phi = Mh(R_field,U,Ud,Comm);

% 2. Conjugate fields
hh=(randn(Nadj,V,D) + i * randn(Nadj,V,D))/sqrt(2);
H=zeros(Nc,Nc,V,D);
for a=1:Nadj
  H = H + tensor(ge(:,:,a),hh(a,:,:));
end


% 3. Initial half step in position U(dt/2) = exp(i  dt/2 H(0)) U(0)
Dt=i * dt * (0.5);
for mu=1:D
  for x=1:V
    U(:,:,mu,x)=expm( Dt * H(:,:,mu,x)) * U(:,:,mu,x);
  end
end
U = reshape(exchange_fields(reshape(U,Nc*Nc*D,EV),Comm,2),Nc,Nc,D,EV);
Ud=Udag(U,Comm);

X=Phi;
for k=0:N_MD-1 %%%%%%%%%%% start Molecular Dynamics loop.
  
  % 4. Compute the Force in the middle points 
  % 4.a first prepare the gauge part
  
  stap=staples_1x1(U,Comm,0);
  dh = zeros(Nc,Nc);
  beta_ov_nc=beta/Nc;
  
  % 4.b then Compute X = inv(MdM) Phi
  MdM = @(p_)Mh(Mh(p_,U,Ud,Comm),U,Ud,Comm); % handel to the fermion matrix squared.
  X = cgs(@MdM,Phi,tol,maxit,[],[],X);  %% Need Parallel solvers!!!
  MhX = Mh(X,U,Ud,Comm);

  U1=U;
  U1d=Ud;
  for mu=1:D
    for x=1:V
      % (gauge force)
      dh = beta_ov_nc * U(:,:,mu,x) * stap(:,:,mu,x);

      % (fermionic force)
      comloc.extind=???; %% SEGNO: still not clear how to restrit the action of Mh (see: ferm_ind_for_force...)
      xup= Comm.neigh(x,mu,UP);
      for a = 1:Nadj
	U1(:,:,mu,x) = uge(:,:,a) * U(:,:,mu,x);
	U1d(:,:,mu,xup) = Ud(:,:,mu,xup) * uge(:,:,a)';
	temp = (conj(X(:,:,comloc.extind)) .* (Mh(MhX,U1,U1d,comloc) - Mh(MhX,U,Ud,comloc))/epsilon);
	dh=dh + ge(:,:,a) * 2.0 * real(temp);
      end
      
      % traceless anti-hermitean part
      dh=0.5(dh - dh');
      dH(:,:,mu,x) = dh - trace(dh)/Nc;
    end
  end
  
  % 4.c Once the force is computed, add it to H(kdt+dt) = H(kdt) - dt F(kdt + dt/2)
  H = H + dH;
  
  if(k<N_MD-1) 
    % 5.a Unless this is the final step compute U at the next intermediate point.
    %     U(kdt + 3/2 dt) = exp(i dt H(kdt +dt)) U(kdt + dt/2)) 
    Dt=i * dt;
    for mu=1:D
      for x=1:V
	U(:,:,mu,x)=expm( Dt * H(:,:,mu,x)) * U(:,:,mu,x);
      end
    end
  else
    % 5.b Otherwise, finish up neatly by computing U(N dt).
    %     U(Ndt) = exp(i dt/2 H(Ndt)) U((N-1/2) dt))
    Dt=i * dt/2;
    for mu=1:D
      for x=1:V
	U(:,:,mu,x)=expm( Dt * H(:,:,mu,x)) * U(:,:,mu,x);
      end
    end
  end % if (k<N_MD-1)
  U = reshape(exchange_fields(reshape(U,Nc*Nc*D,EV),Comm,2),Nc,Nc,D,EV);
  Ud=Udag(U,Comm);

end % End Molecular Dynamics loop     for k=0:N_MD-1
