function [Iout_ferm Iout_gauge] = fermion_ind_for_force_in_U(ix,Iin,Comm,indexord)

%%% SEGNO (tutto da ripensare cosa voglio)
% una prima parte definisce I0_ferm associato a ix=1 (che potra' essere automatizzata in futuro)
% una seconda parte costruisce Iout_ferm(ix,:) per tutti gli altri ix
% Poi stessa cosa per Iout_gauge

if ( Iin == [] ) % initialize the index
  UP=1;
  DN=2;
  EL=Comm.el;
  D=length(EL);
  Ext_Ind=Comm.extind;
  Iout = ix
  ecb=coord_basis(EL,indexord);

  for mu=1:D
    p_up=Neigh(ix,mu,UP);
    p_dn=Neigh(ix,mu,DN);
    Iout=union(Iout,p_up,p_dn);
  end

else

  x=index2coord(ix);
  for iy=I
    iy_shifted=coord2index(mod(index2coord(iy,ecb,indexord)+x,EL),ecb)
    
    Iout= union(Iout,iy_shifted)

  end
  

end