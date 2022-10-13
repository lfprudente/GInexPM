function [W] = spec_proj(Y)

% Compute the complete Spectral Decomposition of Y

[V,D] = eig(full(Y)); 

% Project the eigenvalues in the unit simplex

[sigma,k] = simplex_proj( diag(D) );

% Build the projection

W = V * diag(sigma) * V';

W = ( W + W' ) / 2.0;