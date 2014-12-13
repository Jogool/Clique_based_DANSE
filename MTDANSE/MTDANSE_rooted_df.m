function [node] = MTDANSE_rooted_df(node,ii)
% given a sink (root) node and pre-existing tree, find the data flow away 
% from the root node
nb_nodes = size(node,2);
dim_DANSE = node(1).dimDANSE;

for jj = node(ii).ff_rec
    idx = node(ii).ff_rec;
    idx = idx(find(idx~= jj));
    % if receiving node is in the same clique as the broadcasting node
    if and(node(jj).isclq,eq(node(jj).clq,node(ii).clq))
        idx = setdiff(idx,node(ii).clq_nbrs);
    else 
        idx = unique([idx setdiff(node(ii).clq_nbrs,node(ii).ff_trans)]);
    end
    
    % there are times when the root node does not have any other signals to
    % add, i.e., a line topology when the root is at one of the ends,
    % therefore add this if statement to catch if this happens
    if isempty(idx)
        if isempty(node(ii).ff_trans)
            node(ii).df(jj).zx = node(ii).loc_zx;
            node(ii).df(jj).zn = node(ii).loc_zn;
        else
            z_x_seq = node(node(ii).ff_trans).df(ii).zx;
            z_n_seq = node(node(ii).ff_trans).df(ii).zn;
            
            gkq_coeff = node(ii).gkq(node(ii).ff_trans).coeff;
            
            node(ii).df(jj).zx = node(ii).loc_zx + (gkq_coeff'*z_x_seq')';
            node(ii).df(jj).zn = node(ii).loc_zn + (gkq_coeff'*z_n_seq')';
        end
    else
        if isempty(node(ii).ff_trans)
            z_x_seq = [node(idx).ff_zx];
            z_n_seq = [node(idx).ff_zn];    
        else
            z_x_seq = [node(node(ii).ff_trans).df(ii).zx node(idx).ff_zx];
            z_n_seq = [node(node(ii).ff_trans).df(ii).zn node(idx).ff_zn];
        end
        
    % get gkq coefficeints used for fusion
    gkq_coeff = [node(ii).gkq(node(ii).ff_trans).coeff node(ii).gkq(idx).coeff];
    gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
    gkq_coeff = cat(1,gkq_coeff{:});
    
    node(ii).df(jj).zx = node(ii).loc_zx + (gkq_coeff'*z_x_seq')';
    node(ii).df(jj).zn = node(ii).loc_zn + (gkq_coeff'*z_n_seq')';
    
    end
    
    node = MTDANSE_rooted_df(node,jj);
    
end
