% define GL's, gauge-group, gauge-rep, chem-pot. gauge coupling, mass(es)
% (regularization), (background-metric), 

GL=[10 10 10 16];
Nodes=[2 2 1 1];
Nc=3;
Nf=2;
beta=6.0;
mass=0.01;
N_traj=10;
N_MD=10;
metric=[1 1 1 1];
init_flag=1; %1: start from random (with spread rnd_gauge_sprd); -1 readin file gauge_conf_file.
rnd_gauge_sprd=0.3;
gauge_conf_file=[];
FermMatType=0;
time_dir=4;
D=length(GL);
indexord=fliplr([1:D]); % here I decide which is the slowest and fastest direction in loops.
thick=1; % thickness of the borders (keep it to one except for dramatic reasons).

%%%%%%%%%%%%%%%%%% SET PARAMETERS ONLY ABOVE THIS LINE

% variables
GV=prod(GL);
T=GL(time_dir);
dt=0.5/N_MD;
S=2^(floor(D/2));
Nadj=Nc^2-1;
ge=sun_gen(Nc,2);

% init nearest neighbour (and MPI if necessary)
Comm = init_geometry(L,Nodes,indexord,thick);
V=prod(L);

% handel for the hermitian fermion matrix
Mh = @(p_,U_,Ud_,Comm_)fermion_matrix_herm(p_,U_,Ud_,Comm_,Nc,mass,gi,gv);

% init gauge 
U = init_gauge_fields(Nc,init_flag,rnd_gauge_sprd,gauge_conf_file,Comm);
Ud=Udag(U,Comm);

% evolution algorithm
for j=1:N_traj
 
 U = MolecularDynamics(U,@Mh,beta,ge,Nc,Nf,dt,N_MD,Comm);
 %HybridMonteCarlo
 %Langevin

end % j=1:N_traj
