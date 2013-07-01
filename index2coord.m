function [coord] = index2coord(index,c_basis,indexord)
% [coord] = index2coord(index,c_basis,indexord)
%
% Translates sequential indices (= integers from 1 to V=prod(L))
% into coordinates (=vectors whose i-th component goes from 0 to L(mu)-1).
% indexord is the list of all directions from fastest to slowest 
% and must be the one used to compute c_basis=coord_basis(L,indexord)
% 
% see also init_geometry, coord_basis, coord2index

index=index-1;
for mu=fliplr(indexord)
 coord(mu) = floor(index/c_basis(mu));
 index = mod(index,c_basis(mu));
end
