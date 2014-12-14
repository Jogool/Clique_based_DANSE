function [node,uo,uo_clq] = construct_tree_clq(node)
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
%   uo_clq          -   update order for MTDANSE
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

nb_clqs = max([node.clq]);
non_clq_node_idx = find([node.clq] == 0);
nb_non_clq_nodes = numel(non_clq_node_idx);

A_clq = zeros(nb_clqs+nb_non_clq_nodes,nb_clqs+nb_non_clq_nodes);
D_clq = A_clq;
A_clq(end - nb_non_clq_nodes + 1:end,end - nb_non_clq_nodes + 1:end) = ...
    A(non_clq_node_idx,non_clq_node_idx);
D_clq(end - nb_non_clq_nodes + 1:end,end - nb_non_clq_nodes + 1:end) = ...
    D(non_clq_node_idx,non_clq_node_idx);
% construct small adjacency matrix that represents all nodes in a clique as
% one index in the adjacency matrix
for ii = 1:nb_clqs    
    clq_idx = find([node.clq] == ii);
    idx = setdiff(1:nb_clqs,ii);
    for jj = idx
        non_clq_idx = find([node.clq] == jj);
        if ~isempty(find(A(clq_idx,non_clq_idx)))
            A_clq(ii,jj) = 1;
            A_clq(jj,ii) = 1;
            [C_clq,~] = min(min(D(clq_idx,non_clq_idx)));
            %[C_non_clq,I_non_clq] = min(min(D(clq_idx,non_clq_idx)));
            D_clq(ii,jj) = C_clq;
            D_clq(jj,ii) = C_clq;
        end
    end
    for jj = 1:nb_non_clq_nodes
        if ~isempty(find(A(clq_idx,non_clq_node_idx(jj))))
            A_clq(ii,nb_clqs+jj) = 1;
            A_clq(nb_clqs+jj,ii) = 1;
            [C_clq,~] = min(min(D(clq_idx,non_clq_node_idx(jj))));
            D_clq(ii,nb_clqs+jj) = C_clq;
            D_clq(nb_clqs+jj,ii) = C_clq;
        end
    end
end

A_tril = tril(A_clq);
euc_weight = D_clq(find(A_tril));
[rows,cols] = find(A_tril);

UG = sparse(rows,cols,euc_weight,nb_clqs+nb_non_clq_nodes,nb_clqs+nb_non_clq_nodes);

[ST,~] = graphminspantree(UG);

A_mst_clq = full(ST)+full(ST)';    % adjanceny matrix of tree topology
A_mst_clq(find(A_mst_clq)) = 1;  

A_mst = zeros(nb_nodes);
A_mst(non_clq_node_idx,non_clq_node_idx) = ...
    A_mst_clq(end - nb_non_clq_nodes + 1:end,end - nb_non_clq_nodes + 1:end);
for ii = 1:nb_clqs
    clq_idx = find([node.clq] == ii);
    idx = setdiff(1:nb_clqs,ii);
    % find minimum distance between connected cliques
    for jj = idx
        if A_mst_clq(ii,jj)
            non_clq_idx = find([node.clq] == jj);
            D_temp = D(clq_idx,non_clq_idx);
            [~,I_clq] = min(min(D_temp,[],2));
            [~,I_non_clq] = min(min(D_temp));
            A_mst(clq_idx(I_clq),non_clq_idx(I_non_clq)) = 1;
            A_mst(non_clq_idx(I_non_clq),clq_idx(I_clq)) = 1;
        end
    end
    % find mimumn distance between clique and non-clique node
    for jj = 1:nb_non_clq_nodes
        if A_mst_clq(ii,nb_clqs+jj)
            D_temp = D(clq_idx,non_clq_node_idx(jj));
            [~,I_clq] = min(D_temp);
            A_mst(clq_idx(I_clq),non_clq_node_idx(jj)) = 1;
            A_mst(non_clq_node_idx(jj),clq_idx(I_clq)) = 1;
        end
    end
end
% connect cliques into a tree
for ii = 1:nb_clqs
    clq_idx = find([node.clq] == ii);
    for jj = 1:numel(clq_idx)
        kk = setdiff(clq_idx,clq_idx(jj));
        [~,I_clq] = min(D(clq_idx(jj),kk));
        A_mst(clq_idx(jj),kk(I_clq)) = 1;
        A_mst(kk(I_clq),clq_idx(jj)) = 1;
    end
end

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

% remove redudant updating from clique path : removes all clique only 
% nodes (only let them update once)

uo_clq = uo;
idx = find([node.isclq]);
for ii = 1:numel(idx) 
    if isempty(setdiff(node(idx(ii)).clq_conn,node(idx(ii)).clq_nbrs))
        [~,I] = find(uo_clq == idx(ii));
        uo_clq(I(2:end)) = [];
    else
        [~,I] = find(uo_clq == idx(ii));
        idx_rem = [];
        for jj = 1:numel(I)
            % checK to make sure that node is not at the end or the
            % beginning of the update path
            if and(~eq(I(jj),1),lt(I,numel(uo_clq)))
                start_node = uo_clq(I(jj)-1);
                end_node = uo_clq(I(jj)+1);
                [~,~,ib] = intersect(node(idx(ii)).clq_nbrs,[start_node end_node]);
                if eq(numel(ib),2)
                    idx_rem = [idx_rem I(jj)];
                end
            end
        end
        uo_clq(idx_rem) = [];
    end
end

