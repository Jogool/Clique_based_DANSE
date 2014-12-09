function [node,uo] = construct_tree_clq(node)
% constrcut a MST based on current node positions
%
% The function generates a MST based on current node positions
%
% Syntax:  [node] = network_gen_tree(node)
%
% Inputs:
%   node            -   contains node information in structure format
%
% Outputs:
%   node            -   contains node with updated tree connections 
%                       in structure format
%   A               -   adjacency matrix of tree
%   uo              -   update order for TDANSE
%                   
% Example:
%    [node,~,~] = construct_tree(node)
%
% Other m-files required: path_find

% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% October 2014; Last revision: 04-Dec-2014

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

% truncate connections such that non-clique nodes only communicate with the
% nearest clique node neighbor
A_clq = A;
for ii = 1:max([node.clq])
    idx = find([node.clq] == ii); % nodes in clique
    idx1 = 1:nb_nodes;
    idx1(idx) = [];                         % nodes not in clique
    idx2 = find(sum(A_clq(idx,idx1)) > 1);
    if ~isempty(idx2)
        for jj = 1:numel(idx2)
            [~,I] = min(D(idx,idx1(idx2(jj))));
            A_clq(idx(idx~=idx(I)),idx1(idx2(jj))) = 0;
            A_clq(idx1(idx2(jj)),idx(idx~=idx(I))) = 0;            
        end
    end
end

% find miminum spanning tree based on current connected adjanceny matrix
% Euclidean weights
A_tril = tril(A_clq);
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
end
