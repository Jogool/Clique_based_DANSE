function [node] = TDANSE(node,node_update)
% The function performs the T-DANSE algorithm in a sequential fashion
%
% Syntax:  [node] = TDANSE(node,node_update)
%
% Inputs:
%   node            -   node structure from network_gen_tree.m
%   node_update     -   which node is updating during this iteration
%
% Outputs:
%   node            -   node structure that contains the cost at each
%                       iteration of the T-DANSE algorithm as well as the
%                       optimal filters
%
% Example:
%    [node] = TDANSE(node,n)
%
% Other m-files required: none
% Subfunctions: node
% MAT-files required: none
%

% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Aug. 2014; Last revision: 16-Sep-2014
[node.cost] = deal(0);
nb_nodes = size(node,2);
dim_DANSE = node(1).dimDANSE;
for ii = 1:nb_nodes
    node(ii).loc_zx = (node(ii).loc_filt_coeff'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).loc_filt_coeff'*node(ii).ss_noise')';
end

node = TDANSE_rooted_ff(node,node_update);
node = TDANSE_rooted_df(node,node_update);
node = TDANSE_filt_update(node,node_update);

%  Calculate cost at each node
for ii=1:nb_nodes
    idx = node(ii).ff_rec;
    z_x_seq = [node(node(ii).ff_trans).df(ii).zx node(idx).ff_zx];
    z_n_seq = [node(node(ii).ff_trans).df(ii).zn node(idx).ff_zn];
    
    gkq_coeff = [node(ii).gkq(node(ii).ff_trans).coeff node(ii).gkq(idx).coeff];
    gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
    gkq_coeff = cat(1,gkq_coeff{:});
    
    % cost at node during current iteration
    node(ii).cost(1) = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        [node(ii).loc_filt_coeff' gkq_coeff']*...
        ([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
end
