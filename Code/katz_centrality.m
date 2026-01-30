function kc = katz_centrality(adj_matrix, alpha, beta, tol, max_iter)
% This function computes the Katz centrality for each node in a graph.
% Input:
% adj_matrix - adjacency matrix of the graph
% alpha - attenuation factor (should be less than 1/max(eigenvalue))
% beta - initial centrality weight (default 1)
% tol - tolerance for convergence (default 1e-6)
% max_iter - maximum number of iterations (default 100)
% Output:
% kc - Katz centrality vector

if nargin < 3
    beta = 1; % Default value for beta
end
if nargin < 4
    tol = 1e-6; % Default tolerance
end
if nargin < 5
    max_iter = 100; % Default max iterations
end

% Number of nodes in the graph
num_nodes = size(adj_matrix, 1);

% Initialize centrality vector
kc = ones(num_nodes, 1) * beta;

% Iterate to compute centrality
for iter = 1:max_iter
    kc_new = alpha * (adj_matrix * kc) + beta;
    
    % Check for convergence
    if norm(kc_new - kc, Inf) < tol
        break;
    end
    
    % Update Katz centrality
    kc = kc_new;
end
end
