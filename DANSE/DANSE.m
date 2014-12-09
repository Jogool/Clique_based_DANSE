function [node] = DANSE(node,node_update)

[node.cost_DANSE] = deal(0);
nb_nodes = size(node,2);
dim_DANSE = node(1).dimDANSE;
% generate broadcast signals
for ii = 1:nb_nodes
    idx = find(ii ~= 1:nb_nodes);
    node(ii).loc_zx = (node(ii).loc_filt_coeff'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).loc_filt_coeff'*node(ii).ss_noise')';
end
    % update external filter coefficients at root node
    node = DANSE_filt_update(node,node_update);
    
    %%  Calculate cost at each node
    for ii=1:nb_nodes
        idx = find(ii ~= 1:nb_nodes);
        z_x_seq = [node(idx).loc_zx];
        z_n_seq = [node(idx).loc_zn];
        
        gkq_coeff = [node(ii).gkq.coeff];
        gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
        gkq_coeff = cat(1,gkq_coeff{:});
        temp_filt = [node(ii).loc_filt_coeff' gkq_coeff'];
        
        % cost at node during current iteration
        node(ii).cost(1) = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        temp_filt*([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
    end
end