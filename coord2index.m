function [index] = coord2index(coord,c_basis)
% [index] = coord2index(coord,c_basis)
%
% Translates coordinates (=vectors whose i-th component goes from 0 to L(i)-1).
% into sequential indices (= integers going from 1 to Vol).
%
% see also init_geometry, coord_basis, index2coord, coord2index

index = dot(coord , c_basis) + 1;
