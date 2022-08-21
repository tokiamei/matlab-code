function psiOut = doApplyHam(psiIn,hloc,N,usePBC)
% function psiOut = doApplyHam(psiIn,hloc,N,usePBC)
% ------------------------
% by Glen Evenbly (c) for www.tensors.net, (v1.1) - last modified 06/2020
% 
% Applies local Hamiltonian (given as sum of nearest neighbor terms 'hloc')
% to input state 'psiIn'. Number of lattice sites specified as 'N' while
% 'usePBC' determines whether open or periodic boundaries are used.

d = size(hloc,1);
psiOut = zeros(d^N,1);
for k = 1:N-1
  % apply local Hamiltonian terms to sites [k,k+1]
  axes = {[2],[2]};
  psi_temp = tensordot(reshape(hloc,[d^2,d^2]),...
    reshape(psiIn,[d^(k-1), d^2, d^(N-1-k)]),axes);
  psiOut = psiOut + reshape(permute(psi_temp,[2,1,3]),[d^N,1]);
end

if usePBC
  % apply periodic term
  axes = {[3,4],[3,1]};
  psi_temp = tensordot(reshape(hloc,[d,d,d,d]),...
    reshape(psiIn,[d, d^(N-2), d]),axes);
  psiOut = psiOut + reshape(permute(psi_temp,[2,3,1]),[d^N,1]);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function C = tensordot(A,B,axes)
% C = tensordot(A,B,axes)
%
% Function for taking the produce of two tensors A and B, designed to mimic the
% numpy tensordot. Input `axes = {A_axes,B_axes}` where `A_axes` and
% `B_axes` are vectors describing the indices to be contracted, e.g. set
% axes = {[2],[1]} to contract the 2nd index of A with the 1st index of B.

A_sz = size(A);
A_inds = 1:ndims(A);
A_inds(axes{1}) = [];
A_shape = [prod(A_sz(A_inds)),prod(A_sz(axes{1}))];

B_sz = size(B);
B_inds = 1:ndims(B);
B_inds(axes{2}) = [];
B_shape = [prod(B_sz(axes{2})),prod(B_sz(B_inds))];

C = reshape(reshape(permute(A,[A_inds,axes{1}]),A_shape) *...
  reshape(permute(B,[axes{2},B_inds]),B_shape),[A_sz(A_inds),B_sz(B_inds)]);