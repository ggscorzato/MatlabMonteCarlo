function U = MolecularDynamicsRPhi(U,@Mh,beta,ge,Nc,Nf,Nadj,dt,N_MD,algo,Comm)
% U = MoleculerDynamicsRPhi(U,@Mh,beta,ge,Nc,Nf,Nadj,dt,N_MD,algo,Comm)
%
% Phi- and R-algorithms as described in Ref. Gottlieb et. al. Phys. Rev. D35 (1987) 2531
% R-algorithm follows exactly the description of pag. 2539 of that paper.
% For the Phi-algorithm I actually write a version which is as near as possible to the  R-algorithm
% (in particular, U lives in midpoints and H in the endpoints of the dt intervals).
% Phi-algorithm with N_MD=1 coincides with a Langevin algorithm, although made of 2 half steps.

L=Comm.l;
D=length(L);
V=prod(L);
S=2^(floor(D/2));
EV=prod(Comm.el);
tol=1e-6;
maxit=1e4;

if (algo=='Phi')
 Nf=0; %in Phi-algo, Nf is not a real parameter. The number of fermions is determined by the fermion matrix
end

% (0) "Now let us begin with U(t) and a newly refreshed H(t)
hh=(randn(Nadj,V,D) + i * randn(Nadj,V,D))/sqrt(2);
H=zeros(Nc,Nc,V,D);
for a=1:Nadj
 H = H + tensor(ge(:,:,a),hh(a,:,:));
end

% (1) generate an internediate U
% R-algo: U(t+ (0.5 - Nf/8) dt) = exp(i (0.5 - Nf/8) dt H(t)) U(t)   Eq. (48)
% Phi-algo: U(t+ dt/2) = exp(i  dt/2 H(t)) U(t)
Dt=i * dt * (0.5 - Nf/8);
for mu=1:D
 for x=1:V
  U(:,:,mu,x)=expm( Dt * H(:,:,mu,x)) * U(:,:,mu,x);
 end
end
U = reshape(exchange_fields(reshape(U,Nc*Nc*D,EV),Comm,2),Nc,Nc,D,EV);
Ud=Udag(U,Comm);

% (2-Phi) In Phi-algo the Phi field is generated at the beginning and kept fixed along the MD evolution 
if (algo=='Phi') 
 R_field = (randn(Nc,S,V) + i * randn(Nc,S,V))/sqrt(2); % normalization?
 Phi = Mh(R_field,U,Ud); % Eq. (7)
end

X=Phi;
for n=1:N_MD %%%%%%%%%%% start Molecular Dynamics loop.

 if (algo=='R') % In R-algo the Phi field is generated new at each MD step.
  % (2-R) Generate a gaussian random vector R and an intermediate vector \Phi. 
  % with distribution  P(\Phi)=e^{-\Phi^d (M^dM)^-1 \Phi}  
  R_field = (randn(Nc,S,V) + i * randn(Nc,S,V))/sqrt(2); % normalization?
  Phi = Mh(R_field,U,Ud); % Eq. (49)

  % (3-R) compute U at the midpoint
  % U(t + dt/2) = exp(i (Nf/8) dt H(t)) U(t + (0.5 - Nf/8) dt)   Eq. (50)
  Dt=i * dt * Nf/8;
  for mu=1:D
   for x=1:V
     U(:,:,mu,x)=expm( Dt * H(:,:,mu,x)) * U(:,:,mu,x);
   end
  end
 end
 U = reshape(exchange_fields(reshape(U,Nc*Nc*D,EV),Comm,2),Nc,Nc,D,EV);
 Ud=Udag(U,Comm);

 % (4) Compute X
 % X = inv(MdM) \Phi Eq. (51))
 MdM = @(p_)Mh(Mh(p_,U,Ud),U,Ud); % handel to the fermion matrix squared.
 X= cgs(@MdM,Phi,tol,maxit);  %% Need Parallel solvers!!!
			      %% start from previous solution as a guess if Phi algo is used. !!! 
			      %%  X = CGS(@MdM,Phi,tol,maxit,[],[],X)

 % (5) Use Eq. (46) and (47) to compute dH/dt
 stap=staples_1x1(U,Comm,0);
 dh = zeros(Nc,Nc);
 beta_ov_nc=beta/Nc;
 for mu=1:D
  for x=1:V
   % (5-gauge)
   dh = beta_ov_nc * U(:,:,mu,x) * stap(:,:,mu,x);
   % (5-ferm)   
   %%%%%%%%%% SEGNO
   dh=dh + ??
   % traceless anti-hermitean part
   dh=0.5(dh - dh');
   dH(:,:,mu,x) = dh - trace(dh)/Nc;
  end
 end

 % (6) Compute H(t+dt)
 H = H + dH; % Eq. (52)

 % (7) 
 if(n<N_MD) % Unless this is the last time step ...
	    % ... compute U at the next intermediate point.
	    % R-algo: U(t + dt + (0.5-Nf/8) dt) = exp(i (1 - Nf/8) dt H(t+dt)) U(t+dt/2))   Eq. (53)
	    % Phi-algo: U(t + dt/2 + dt) = exp(i dt H(t+dt)) U(t+dt/2)) 
  Dt=i * dt * (1.0 - Nf/8);
  for mu=1:D
   for x=1:V
    U(:,:,mu,x)=expm( Dt * H(:,:,mu,x)) * U(:,:,mu,x);
   end
  end
  U = reshape(exchange_fields(reshape(U,Nc*Nc*D,EV),Comm,2),Nc,Nc,D,EV);
  Ud=Udag(U,Comm);
 else
  % Otherwise finish up neatly by computing U(t+N dt).
  % U(t + N dt) = exp(i dt/2 H(t+dt)) U(t + (N-1/2) dt))   Eq. (54)
  Dt=i * dt/2;
  for mu=1:D
   for x=1:V
    U(:,:,mu,x)=expm( Dt * H(:,:,mu,x)) * U(:,:,mu,x);
   end
  end
  U = reshape(exchange_fields(reshape(U,Nc*Nc*D,EV),Comm,2),Nc,Nc,D,EV);
  Ud=Udag(U,Comm);
 end % if (n<N_MD)
end %%%%%%%%%%%%% End Molecular Dynamics loop     for n=1:N_MD
