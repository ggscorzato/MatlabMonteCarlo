GL=[4 4 4 8]; % sizes of the global physical lattice
Nodes=[1 1 1 1]; % number of processors in which each direction is subdivided.
Nc=3; % su(Nc) gauge group
Nf=2; % number of flavours
beta=6.0;
mass=0.01;
N_traj=10;
N_MD=10;
metric=[1 1 1 1];
init_flag=1; %1: start from random (with spread rnd_gauge_sprd); -1 readin file gauge_conf_file.
rnd_gauge_sprd=0.3;
gauge_conf_file=[];
time_dir=4;
D=length(GL);
indexord=fliplr([1:D]); % here I decide which is the slowest and fastest direction in loops.

%%%%%%%%%%%%%%%%%% SET PARAMETERS ONLY ABOVE THIS LINE

% variables
GV=prod(GL);
T=GL(time_dir);
dt=0.5/N_MD;
S=2^(floor(D/2));
Nadj=Nc^2-1; 
ge=sun_gen(Nc,2);

% init nearest neighbour (and MPI if necessary)
comm = init_geometry(L,Nodes,indexord,1);

% init gauge 
U = init_gauge_fields(Nc,init_flag,rnd_gauge_sprd,gauge_conf_file,comm);

% communications
