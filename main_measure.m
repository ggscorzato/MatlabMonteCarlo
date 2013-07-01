% define L's, gauge-group, gauge-rep, chem-pot. gauge coupling, mass(es)
% (regularization), (background-metric), 
GL=[4 4 4 4];
Nodes=[2 2 1 1];
Nc=3;
mass=0.1;
N_meas=1;
alpha=0.01;
D=length(L);
V=prod(L);
S=2^(floor(D/2));
indexord=fliplr([1:D]); % here I decide which is the slowest and fastest direction in loops.

% init nearest neighbours  (geometry)
[Neigh Ext_Ind L EL] = init_geometry(GL,Nodes,indexord);

% build SU(Nc) generators.
ge=sun_gen(Nc,2);  %%%

% build D dimensional gamma matrices.
[gi gv ] = cliff_gen(0,D);  %%

% start loop over configurations to be measured
for j=1:N_meas

% read gauge configuration U(t0)
gauge_init_flag=-1; % random
[U Ud]= init_gauge_fields(ge,L,gauge_init_flag,alpha,'gauge.mat');

% gauge measurements (Wloops, Polyakov)

% test
ss=[1, 1, 57]
spin_init_flag=2;
pp = init_spinor_field(Nc,L,spin_init_flag,ss,[],[]);
for j=1:20
  pp= fermion_matrix(pp,L,Nc,mass,2);
end
for j=1:20
  pp= fermion_matrix(pp,L,Nc,mass,0);
end
for j=1:20
  pp= fermion_matrix(pp,L,Nc,mass,1);
end
save pp

return

% compute propagator 
tol=1e-6;
maxit=200;
modelfun=@(p)reshape(fermion_matrix(reshape(p,Nc,S,V),L,Nc,mass,TYPE),Nc*S*V,1);


%ss=[1, 1, ceil(rand*V)]
ss=[1, 1, 57]
spin_init_flag=2;
for is=1:1 % S
  for ic=1:1 % Nc
    ss(1:2) = [ic,is];
    source = reshape(init_spinor_field(Nc,L,spin_init_flag,ss,[],[]),Nc*S*V,1);
    temp = bicgstab(modelfun,source,tol,maxit,[],[],source);
    prop(ic,is,:,:,:)=reshape(temp,Nc,S,V);
  end
end


% compute correlators

% save results
save gauge.mat U 
save prop.mat prop
end

