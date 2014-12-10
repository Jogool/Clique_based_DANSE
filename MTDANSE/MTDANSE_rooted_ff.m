function [node] = MTDANSE_rooted_ff(node,root)
% given a sink (root) node and pre-existing tree, find the data flow toward the
% sink (root) node
nb_nodes = size(node,2);
dim_DANSE = node(1).dimDANSE;

[node.ff_trans] = deal([]);   % node k transmits to this node during the ff (should always be a single node)
[node.ff_rec] = deal([]);     % node k receives these signals during the ff
[node.ff_update] = deal(0);   % flag if node has transmitted its ff-signal (0-no, 1-yes)

% if node has only clique neighbors in same clique, can generate its
% ff-sginal
clq_idx = find([node.isclq]);
for ii = 1:numel(clq_idx)
    if isempty(setdiff(node(clq_idx(ii)).clq_conn,node(clq_idx(ii)).clq_nbrs));
        % update ff signal
        node(ii).ff_zx = node(ii).loc_zx;
        node(ii).ff_zn = node(ii).loc_zn;
        node(ii).ff_update = 1;
    end
end

% if node only has one connection and is not the root node, then it can
% immediately transmit its fusion flow signal
for ii = find(cellfun(@(x) numel(x), {node.clq_conn}) == 1)
    if and(~eq(ii,root),~node(ii).isclq)
        % which node the current node transmits to during the ff
        node(ii).ff_trans = node(ii).clq_conn;
        % the node that receives the ff from the node
        node(node(ii).clq_conn).ff_rec = ...
            sort([node(node(ii).clq_conn).ff_rec ii]);
        
        % update ff signal
        node(ii).ff_zx = node(ii).loc_zx;
        node(ii).ff_zn = node(ii).loc_zn;
        node(ii).ff_update = 1; 
    end
end

% number of nodes who have performed ff update
node_ff_update = numel(find(cat(1,node.ff_update)));   

while lt(node_ff_update,nb_nodes-1) % can skip the root node, hence - 1
    % find all nodes who have not performed a fusion flow update
    ff_idx = find(~cat(1,node.ff_update));  
    % remove the root node as it will never generate a ff signal
    ff_idx(find(ff_idx == root)) = [];
    for ii = ff_idx'

        % find neighbors of node who have already transmitted a ff signal
        idx = find([node(node(ii).clq_conn).ff_update]);
        nbrs_updated = node(ii).clq_conn(idx);
        nb_nbrs_updated = numel(nbrs_updated);
        
        % find neighbors who have not performed a fusion flow update
        idx = find(~[node(node(ii).clq_conn).ff_update]);
        non_update_neighbors = node(ii).clq_conn(idx);
        
        % if all neighbors except 1 have performed a fusion flow
        % update, the node can generate its fusion flow signal
        if eq(nb_nbrs_updated,numel(node(ii).clq_conn)-1)    
            
            % if we are broadcasting to a node in the same clique, do not
            % include other clique signals
            if and(node(ii).isclq,eq(node(ii).clq,node(non_update_neighbors).clq))
                idx = setdiff(nbrs_updated,node(ii).clq_nbrs);
                
                % gather all updated non-clique neighbor signals
                z_x_seq = [node(idx).ff_zx];
                z_n_seq = [node(idx).ff_zn];
                gkq_coeff = [node(ii).gkq(idx).coeff];
                gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
                gkq_coeff = cat(1,gkq_coeff{:});
                
                 % add local transmitted signals to node's ff signal

                node(ii).ff_zx = node(ii).loc_zx + (gkq_coeff'*z_x_seq')';
                node(ii).ff_zn = node(ii).loc_zn + (gkq_coeff'*z_n_seq')';

            else
                % place node in sorted list of received ff signals of
                % non-updated node
                idx = sort([node(non_update_neighbors).ff_rec ii]);
                node(non_update_neighbors).ff_rec = idx;
                
                node(ii).ff_trans = non_update_neighbors;
                
                % gather all updated neighbor signals
                z_x_seq = [node(nbrs_updated).ff_zx];
                z_n_seq = [node(nbrs_updated).ff_zn];
                
                gkq_coeff = [node(ii).gkq(nbrs_updated).coeff];
                gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
                gkq_coeff = cat(1,gkq_coeff{:});
                
                % add local transmitted signals to node's ff signal
                node(ii).ff_zx = node(ii).loc_zx + (gkq_coeff'*z_x_seq')';
                node(ii).ff_zn = node(ii).loc_zn + (gkq_coeff'*z_n_seq')';
            end
            node(ii).ff_update = 1;
        end
    end
    node_ff_update = numel(find(cat(1,node.ff_update)));
end