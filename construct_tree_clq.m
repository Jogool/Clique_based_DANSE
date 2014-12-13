function [node,uo] = construct_tree_clq(node)
% constrcut a MST based on current node positions
%
% The function generates a MST based on current node positions
%
% Syntax:  [node] = construct_tree_clq(node)
%
% Inputs:
%   node            -   contains node information in structure format
%
% Outputs:
%   node            -   contains node with updated tree connections 
%                       in structure format
%   uo              -   update order for TDANSE
%                   
% Example:
%    [node,~] = construct_tree_clq(node)
%
% Other m-files required: path_find

% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% October 2014; Last revision: 10-Dec-2014

nb_nodes = size(node,2);
r = 0;                              % initial broadcast radius of each node
alg_conn = 0;                       % algebratic connectivity
while lt(alg_conn,1e-10)

    A = zeros(nb_nodes);            % adjacency matrix based on current connectivity      
    D = zeros(nb_nodes);            % matrix that holds the distances between nodes
    for ii = 1:nb_nodes
        index = find(ii ~= 1:nb_nodes);
        for jj = 1:nb_nodes-1
            D(ii,index(jj)) = norm(node(ii).pos-node(index(jj)).pos);
            if lt(D(ii,index(jj)),r)    % if the distance is less than the sensing radius, assume nodes can communicate
                A(ii,index(jj)) = 1;
            end
        end
    end

    Deg_mat = diag(sum(A,2));   % diagnomal matrix is equal to connections
    L = Deg_mat-A;              % Laplacian matrix 
    lambda = eig(L);            % eigenvalues of Laplacian matrix
    [~,I] = sort(lambda);
    alg_conn = lambda(I(2));    % algebratic connectivity is equal to the second largest eigenvalue of Laplacian matrix
    r = r + 0.01;               % increase sensing radius 
end

% reduce node communication range to minimal distance for connectivity
max_range = zeros(nb_nodes,1);
for ii = 1:nb_nodes
    index = find(A(ii,:));
    for jj = 1:length(index)
        max_range(ii) = max(max_range(ii),  D(ii,index(jj)));
    end
end

%find all posssible connections (ad-hoc)
for ii = 1:nb_nodes
    node(ii).conn = find(A(ii,:));
end

% find which clique nodes are closest to non-clique nodes, use these as
% neighbors
non_clq_nodes_idx = find([node.clq] == 0);
nb_clqs = max([node.clq]);
AA = zeros(nb_nodes);
AA(non_clq_nodes_idx,non_clq_nodes_idx) = A(non_clq_nodes_idx,non_clq_nodes_idx);

for ii = 1:nb_clqs
    % this will connect to a whole clique still, perhaps this is the
    % problem?
    clq_idx = find([node.clq] == ii);
    AA(clq_idx,clq_idx) = A(clq_idx,clq_idx);
    non_clq_idx = setdiff(1:nb_nodes,clq_idx);
    [~,non_clq_conn] = find(A(clq_idx,non_clq_idx));
    non_clq_idx = non_clq_idx(unique(non_clq_conn));
    for jj = 1:numel(non_clq_idx)
        D_temp = D(clq_idx,non_clq_idx(jj));
        D_temp(find(D_temp == 0)) = inf;
        [~,I] = min(D_temp);
        AA(clq_idx(I),non_clq_idx(jj)) = 1;
        AA(non_clq_idx(jj),clq_idx(I)) = 1;
    end
end

A_tril = tril(AA);
euc_weight = D(find(A_tril));
[rows,cols] = find(A_tril);

UG = sparse(rows,cols,euc_weight,nb_nodes,nb_nodes);

[ST,~] = graphminspantree(UG);

A_mst = full(ST)+full(ST)';    % adjanceny matrix of tree topology
A_mst(find(A_mst)) = 1;   


% find all connections for tree
for ii = 1:nb_nodes
    node(ii).tree_conn = find(A_mst(:,ii));
end

% add clique connections to tree connections for clique nodes
for ii = 1:nb_nodes
    if node(ii).isclq
        idx = find([node.clq] == node(ii).clq);
        clq_nbrs = idx(find(idx ~= ii));
        node(ii).clq_nbrs = clq_nbrs;
        node(ii).clq_conn = unique([clq_nbrs node(ii).tree_conn'])';
    else
        node(ii).clq_conn = node(ii).tree_conn;
    end
end

[u1,v1]=eig(A_mst);
[~,lambdaind]=sort(abs(diag(v1)));
EigvCentrality=[[1:size(A_mst,1)]' abs(u1(:,lambdaind(end)))];

[~,root_node] = max(EigvCentrality(:,2));   % pick root node based on eigenvalue centrality

% find updating order based on root node and 
uo = root_node;
uo = path_find(uo,A_mst,root_node);

