function hc = harmonic_centrality(adj_matrix)
% This function computes the harmonic centrality for each node in a graph.
% Input: adj_matrix - adjacency matrix of the graph
% Output: hc - harmonic centrality vector

% Number of nodes in the graph
num_nodes = size(adj_matrix, 1);

% Initialize the harmonic centrality vector
hc = zeros(num_nodes, 1);

% Calculate the shortest path distances between all pairs of nodes
dist_matrix = graphallshortestpaths(sparse(adj_matrix));

% Compute harmonic centrality for each node
for v = 1:num_nodes
    % Avoid division by zero for self-loops, hence diag(dist_matrix) should be Inf
    dist = dist_matrix(v, :);
    dist(v) = Inf;  % Set distance to itself as infinity to avoid self-loops
    hc(v) = sum(1 ./ dist(dist ~= Inf));  % Sum of reciprocals of distances
end
end
