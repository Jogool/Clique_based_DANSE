%% update filter coefficents of node
function [node] = DANSE_filt_update(node,ii)
nb_nodes = size(node,2);
dim_DANSE = node(1).dimDANSE;

idx = find(ii ~= 1:nb_nodes);
z_x_seq = [node(idx).loc_zx];
z_n_seq = [node(idx).loc_zn];


Rxx = [node(ii).ss_clean z_x_seq]'*[node(ii).ss_clean z_x_seq];
Rnn = [node(ii).ss_noise z_n_seq]'*[node(ii).ss_noise z_n_seq];

w_temp  = (Rnn+Rxx) \ Rxx(:,1:dim_DANSE);       % update node-specific filter
node(ii).loc_filt_coeff = w_temp(1:node(ii).sensors,:);
gkq_seq_temp = w_temp(node(ii).sensors+1:end,:);

for jj = 1:numel(idx);
    node(ii).gkq(idx(jj)).coeff =  gkq_seq_temp((jj-1)*dim_DANSE+1:jj*dim_DANSE,:);
end