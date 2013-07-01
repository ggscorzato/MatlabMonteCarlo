function [c_basis] = coord_basis(L,indexord)
% [c_basis] = coord_basis(L,indexord)
% 
% this subroutine constructs a vector that can be used to translate indices (= numbers going from 1 to V=prod(L))
% into coordinates (=vectors whose i-th component goes from 0 to L(i)-1) and vice-versa.
% It is used by coord2index() and by index2coord(), but can be used by other routines.
% indexord is any permutation of [1:D]. It is the only thing that must be changed in the code in order to decide 
% which coordinates run slow or fast. The fastest direction is indexord(1) and the slowest indexord(D).
%
% With D=4 and default indexord (= [4 3 2 1]), we have:
% c_basis = [L(2)*L(3)*L(4), L(3)*L(4), L(4), 1]
% ind = x(4) + x(3)*L(4) + x(2)*L(3)*L(4) + x(1)*L(2)*L(3)*L(4) + 1
% 
% see also init_geometry, index2coord, coord2index

D=length(L);

c_basis(indexord(1))=1;
for mu=[2:D]
  c_basis(indexord(mu))=c_basis(indexord(mu-1))*L(indexord(mu-1));
end
