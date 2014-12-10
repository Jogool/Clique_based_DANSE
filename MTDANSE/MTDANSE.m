function [node] = MTDANSE(node,node_update)
% The function performs the DANSE algorithm in a mixed topology
%
% Syntax:  [node] = MTDANSE(node,node_update)
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
%    [node] = MTDANSE(node,n)
%
% Other m-files required: MTDANSE_filt_update, MTDANSE_rooted_df,
% MTDANSE_rooted
%

% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% Dec. 2014; Last revision: 10-Dec-2014
[node.cost] = deal(0);
nb_nodes = size(node,2);
dim_DANSE = node(1).dimDANSE;
for ii = 1:nb_nodes
    node(ii).loc_zx = (node(ii).loc_filt_coeff'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).loc_filt_coeff'*node(ii).ss_noise')';
end

node = MTDANSE_rooted_ff(node,node_update);
node = MTDANSE_rooted_df(node,node_update);
node = MTDANSE_filt_update(node,node_update);



% intra-clqiue signals denoted as clq_zx, clq_zn
% inter-clique and tree signals given as ff_zx, ff_zn, df_zx, df_zn

for ii = 1:nb_nodes
    node(ii).loc_zx = (node(ii).loc_filt_coeff'*node(ii).ss_clean')';
    node(ii).loc_zn = (node(ii).loc_filt_coeff'*node(ii).ss_noise')';
end

%% fusion flow 
for ii  = dd_flow
    
    if node(ii).isClique
        if isempty(node(ii).ff_tran) %%% Intra-clique broadcast signal
            idx = node(ii).ff_rec;
            % if node does not receive any ff signals then it can fire
            % using only its local signals
            if isempty(idx)
               node(ii).clq_zx = node(ii).loc_zx;
               node(ii).clq_zn = node(ii).loc_zn;
            
            % if node does receive ff signals than it  combines these with
            % its local signals to make clique signal
            else       
                z_x_seq = [node(idx).ff_zx];
                z_n_seq = [node(idx).ff_zn];
                
                gkq_coeff = [node(ii).gkq(idx).coeff];
                gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
                gkq_coeff = cat(1,gkq_coeff{:});
                
                node(ii).clq_zx = node(ii).loc_zx + (gkq_coeff'*z_x_seq')';
                node(ii).clq_zn = node(ii).loc_zn + (gkq_coeff'*z_n_seq')';
            end
        else % Inter-clique broadcast signal
            
            idx = node(ii).ff_rec;
            if isempty(idx)
                clq_nbrs = find([node.clq_idx] == node(ii).clq_idx);
                clq_nbrs(find(clq_nbrs == ii)) = [];
                
                z_x_seq = [node(clq_nbrs).clq_zx];
                z_n_seq = [node(clq_nbrs).clq_zn];
                
                gkq_coeff = [node(ii).gkq(clq_nbrs).coeff];
                gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
                gkq_coeff = cat(1,gkq_coeff{:});
            else
                clq_nbrs = find([node.clq_idx] == node(ii).clq_idx);
                clq_nbrs(find(clq_nbrs == ii)) = [];
                
                z_x_seq = [node(idx).ff_zx node(clq_nbrs).clq_zx];
                z_n_seq = [node(idx).ff_zn node(clq_nbrs).clq_zn];
                
                gkq_coeff = [node(ii).gkq(idx).coeff node(ii).gkq(clq_nbrs).coeff];
                gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
                gkq_coeff = cat(1,gkq_coeff{:});
                
            end
            node(ii).ff_zx = node(ii).loc_zx + (gkq_coeff'*z_x_seq')';
            node(ii).ff_zn = node(ii).loc_zn + (gkq_coeff'*z_n_seq')';
            
        end
    else % if node is not a clique node, then it uses the T-DANSE rules
        idx = node(ii).ff_rec;
        % if node does not receive any ff signals then it is a leaf node
        if isempty(idx)
            node(ii).ff_zx = node(ii).loc_zx;
            node(ii).ff_zn = node(ii).loc_zn;
        % if node is not a leaf node, it must combine all of its neighbors
        % signals with its local signals to generate ff signal
        else
            z_x_seq = [node(idx).ff_zx];
            z_n_seq = [node(idx).ff_zn];
            
            gkq_coeff = [node(ii).gkq(idx).coeff];
            gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
            gkq_coeff = cat(1,gkq_coeff{:});
            
            node(ii).ff_zx = node(ii).loc_zx + (gkq_coeff'*z_x_seq')';
            node(ii).ff_zn = node(ii).loc_zn + (gkq_coeff'*z_n_seq')';
        end
    end
end

%% diffusion flow, opposite order of fusion flow
for ii = fliplr(dd_flow)
    
    if node(ii).isClique
        clq_nbrs = find([node.clq_idx] == node(ii).clq_idx);
        clq_nbrs(find(clq_nbrs == ii)) = [];
        
        if and(~isempty(node(ii).ff_tran),isempty(node(ii).ff_rec)) %%% Intra-clique broadcast signal
            idx = node(ii).ff_tran;
            node(ii).clq_zx = node(ii).loc_zx + (node(ii).gkq(idx).coeff'*node(idx).df(ii).zx')';
            node(ii).clq_zn = node(ii).loc_zn + (node(ii).gkq(idx).coeff'*node(idx).df(ii).zn')';
        end
        
        if and(isempty(node(ii).ff_tran),~isempty(node(ii).ff_rec))  %%% inter-clique signals
            z_x_seq = [node(clq_nbrs).clq_zx];
            z_n_seq = [node(clq_nbrs).clq_zn];
            
            gkq_coeff = [node(ii).gkq(clq_nbrs).coeff];
            gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
            gkq_coeff = cat(1,gkq_coeff{:});
            for jj = node(ii).ff_rec
                node(ii).df(jj).zx = node(ii).clq_zx - (node(ii).gkq(jj).coeff'*node(jj).ff_zx')' + (gkq_coeff'*z_x_seq')';
                node(ii).df(jj).zn = node(ii).clq_zn - (node(ii).gkq(jj).coeff'*node(jj).ff_zn')' + (gkq_coeff'*z_n_seq')';
            end
        end
        
        if and(~isempty(node(ii).ff_tran),~isempty(node(ii).ff_rec))  %%% inter, intra and signals
            ff_idx = node(ii).ff_rec;
            df_idx = node(ii).ff_tran;
            
            z_x_seq = [node(df_idx).df(ii).zx node(ff_idx).ff_zx];
            z_n_seq = [node(df_idx).df(ii).zn node(ff_idx).ff_zn];
            
            gkq_coeff = [node(ii).gkq(df_idx).coeff node(ii).gkq(ff_idx).coeff];
            gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
            gkq_coeff = cat(1,gkq_coeff{:});
            
            node(ii).clq_zx = node(ii).loc_zx + (gkq_coeff'*z_x_seq')';
            node(ii).clq_zn = node(ii).loc_zn + (gkq_coeff'*z_x_seq')';
            
            z_x_seq = [node(clq_nbrs).clq_zx];
            z_n_seq = [node(clq_nbrs).clq_zn];
            
            gkq_coeff = [node(ii).gkq(clq_nbrs).coeff];
            gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
            gkq_coeff = cat(1,gkq_coeff{:});
            
            for jj = node(ii).ff_rec
                node(ii).df(jj).zx = node(ii).clq_zx - (node(ii).gkq(jj).coeff'*node(jj).ff_zx')' + (gkq_coeff'*z_x_seq')';
                node(ii).df(jj).zn = node(ii).clq_zn - (node(ii).gkq(jj).coeff'*node(jj).ff_zn')' + (gkq_coeff'*z_n_seq')';
            end
        end
    else
        for jj = node(ii).ff_rec
            idx = node(ii).ff_tran;
            if isempty(idx)
                % root node broadcast all signals minus FF signal from
                % neighbor
                node(ii).df(jj).zx = node(ii).ff_zx - (node(ii).gkq(jj).coeff'*node(jj).ff_zx')';
                node(ii).df(jj).zn = node(ii).ff_zn - (node(ii).gkq(jj).coeff'*node(jj).ff_zn')';
            else
                % TFC signal = nodes local signal + scaled DF signal
                node(ii).df(jj).zx = node(ii).loc_zx + (node(ii).gkq(idx).coeff'*node(idx).df(ii).zx')';% - (node(ii).gkq(jj).coeff'*node(jj).ff_zx')';
                node(ii).df(jj).zn = node(ii).loc_zn + (node(ii).gkq(idx).coeff'*node(idx).df(ii).zn')';% - (node(ii).gkq(jj).coeff'*node(jj).ff_zn')';
            end
        end
    end
end

% update node-specific filter coefficients at updated node
node = node_filt_update(node,node_update);

% calculate cost at each node
for ii=1:nb_nodes
    if node(ii).isClique
        clq_nbrs = find([node.clq_idx] == node(ii).clq_idx);
        clq_nbrs(find(clq_nbrs == ii)) = [];
        z_x_seq = [node(node(ii).ff_tran).df(ii).zx node(node(ii).ff_rec).ff_zx node(clq_nbrs).clq_zx];
        z_n_seq = [node(node(ii).ff_tran).df(ii).zn node(node(ii).ff_rec).ff_zn node(clq_nbrs).clq_zn];
        
        gkq_coeff = [node(ii).gkq(node(ii).ff_tran).coeff node(ii).gkq(node(ii).ff_rec).coeff node(ii).gkq(clq_nbrs).coeff];
        gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
        gkq_coeff = cat(1,gkq_coeff{:});
    else
        z_x_seq = [node(node(ii).ff_tran).df(ii).zx node(node(ii).ff_rec).ff_zx];
        z_n_seq = [node(node(ii).ff_tran).df(ii).zn node(node(ii).ff_rec).ff_zn];
        
        gkq_coeff = [node(ii).gkq(node(ii).ff_tran).coeff node(ii).gkq(node(ii).ff_rec).coeff];
        gkq_coeff = mat2cell(gkq_coeff, size(gkq_coeff,1), dim_DANSE*ones(1,size(gkq_coeff,2)/dim_DANSE));
        gkq_coeff = cat(1,gkq_coeff{:});
    end

    
    % cost at node during current iteration
    node(ii).cost_mixed_DANSE(1) = norm(node(ii).ss_clean(:,1:dim_DANSE)' - ...
        [node(ii).loc_filt_coeff' gkq_coeff']*...
        ([node(ii).ss_clean z_x_seq]+[node(ii).ss_noise z_n_seq])')^2;
end
end
%% node update
function [node] = node_filt_update(node,ii)
global nb_nodes dim_DANSE

clq_nbrs = [];
if node(ii).isClique
    clq_nbrs = find([node.clq_idx] == node(ii).clq_idx);
    clq_nbrs(find(clq_nbrs == ii)) = [];
    z_x_seq = [node(node(ii).ff_tran).df(ii).zx node(node(ii).ff_rec).ff_zx node(clq_nbrs).clq_zx];
    z_n_seq = [node(node(ii).ff_tran).df(ii).zn node(node(ii).ff_rec).ff_zn node(clq_nbrs).clq_zn];
else
    z_x_seq = [node(node(ii).ff_tran).df(ii).zx node(node(ii).ff_rec).ff_zx];
    z_n_seq = [node(node(ii).ff_tran).df(ii).zn node(node(ii).ff_rec).ff_zn];
end
Rxx = [node(ii).ss_clean z_x_seq]'*[node(ii).ss_clean z_x_seq];
Rnn = [node(ii).ss_noise z_n_seq]'*[node(ii).ss_noise z_n_seq];
w_temp  = pinv(Rnn+Rxx) * Rxx(:,1:dim_DANSE);       % update node-specific filter
node(ii).loc_filt_coeff = w_temp(1:node(ii).sensors,:);
gkq_seq_temp = w_temp(node(ii).sensors+1:end,:);

idx = [node(ii).ff_tran node(ii).ff_rec clq_nbrs];

for jj = 1:numel(idx);
    node(ii).gkq(idx(jj)).coeff =  gkq_seq_temp((jj-1)*dim_DANSE+1:jj*dim_DANSE,:);
end

end
%% have each node find its place in the topology
function [node,ff_update] = node_order(node)
[node.ff_update] = deal(0);
[node.clq_update] = deal(0);
ff_update = [];

% find all nodes which only have clique neighbors, they can update automatically 
idx = find([node.isClique]);
for ii = idx
    clq_nbrs = find([node.clq_idx] == node(ii).clq_idx);
    clq_nbrs(find(clq_nbrs == ii)) = [];
    nb_clq_nbrs = numel(clq_nbrs);
    
    if eq(nb_clq_nbrs,numel(node(ii).clq_connections))
        node(ii).ff_rec = [];
        node(ii).clq_update = 1;
        ff_update = [ff_update ii];
    end
end

% find all leaf nodes, they can update automatically
idx = find(~[node.isClique]);
for ii = idx
    if eq(numel(node(ii).clq_connections),1)
        node(ii).ff_tran = node(ii).clq_connections;
        node(ii).ff_update = 1;
        ff_update = [ff_update ii];
        node(ii).ff_rec = [];       
    end
end
idx = find(~or([node.clq_update],[node.ff_update]));
teller = 1;

% note that when idx is 1 then that node is the root node, if idx is empty,
% the fusion flow ends at a clique with more that one node, i.e., the
% clique can be abstracted as the root node of the system
while gt(numel(idx),1) 
    %or(~eq(numel(idx),1),isempty(idx)) 
    for ii = idx
        if node(ii).isClique
            % find all clique neighbors who only have neighbors in the clique
            %(they have already update their intra-clique signals)
            clq_nbrs = find([node.clq_idx] == node(ii).clq_idx);
            clq_nbrs(find(clq_nbrs == ii)) = [];
            clq_nbrs_updated = clq_nbrs(find([node(clq_nbrs).clq_update]));
            nb_clq_nbrs_updated = numel(clq_nbrs_updated);
            
            % find all non clique-neighbors who have trasmitted their ff-signal
            non_clq_nbrs = setdiff(node(ii).clq_connections,clq_nbrs);
            non_clq_nbrs_updated = non_clq_nbrs(find([node(non_clq_nbrs).ff_update]))';
            nb_non_clq_nbrs_updated = numel(non_clq_nbrs_updated);
           
            nb_nbrs_updated = nb_clq_nbrs_updated + nb_non_clq_nbrs_updated;
            
            % find all neighbors that have not updated
            non_update_neighbors = setdiff(node(ii).clq_connections,[clq_nbrs_updated non_clq_nbrs_updated]);
            nb_non_update_neighbors = numel(non_update_neighbors); 
            
            % check to see if all neighbors that have not updated belong to
            % the same clique as the current node, if they are then this
            % violates looking for only one neighbor that has not updated
            % and the node can essentially flood the clique with its
            % clq_signal
            clq_flag = eq(numel(find([node(non_update_neighbors).clq_idx] == node(ii).clq_idx)),nb_non_update_neighbors);
            if clq_flag
                node(ii).ff_rec = [non_clq_nbrs];
                node(ii).ff_rec = (node(ii).ff_rec(:))';
                ff_update = [ff_update ii];
                node(ii).clq_update = 1;
            end
                
            if and(eq(nb_non_update_neighbors,1),~clq_flag)
                % has to be in same clique, otherwise its an ff signal
                if and(node(non_update_neighbors).isClique,eq(node(non_update_neighbors).clq_idx,node(ii).clq_idx))
                    node(ii).ff_rec = [non_clq_nbrs];

                    ff_update = [ff_update ii];
                    node(ii).clq_update = 1;
                else
                    %node(non_update_neighbors).ff_signals = sort([node(non_update_neighbors).ff_signals ii]);
                    node(ii).ff_tran = non_update_neighbors;
                    % since the clique node is broadcasting to a node not
                    % in the same clique, the ff_signal is composed for not
                    % only the clique nodes but all non clique nodes
                    node(ii).ff_rec = [non_clq_nbrs_updated];
                    node(ii).ff_rec = (node(ii).ff_rec(:))';
                    ff_update = [ff_update ii];
                    node(ii).ff_update = 1;
                end
            end
            
        else
            nb_nbrs_updated = sum([node(node(ii).clq_connections).clq_update])+sum([node(node(ii).clq_connections).ff_update]);
            if eq(nb_nbrs_updated,numel(node(ii).clq_connections)-1)
                idx_tran = find(~or([node((node(ii).clq_connections)).clq_update],[node((node(ii).clq_connections)).ff_update]));
                idx_rec = find(or([node((node(ii).clq_connections)).clq_update],[node((node(ii).clq_connections)).ff_update]));
                ff_update = [ff_update ii];
                node(ii).ff_tran = node(ii).clq_connections(idx_tran);
                node(ii).ff_rec = [node(ii).clq_connections(idx_rec)];
                
                node(ii).ff_tran = (node(ii).ff_tran(:))';
                node(ii).ff_rec = (node(ii).ff_rec(:))';

                node(ii).ff_update = 1;
                
            end
        end
    end
    idx = find(~or([node.clq_update],[node.ff_update]));
    teller = teller + 1;
    if gt(teller,50)
        disp('')
    end 
end
if eq(numel(idx),1) % root node
    if node(idx).isClique
        clq_nbrs = find([node.clq_idx] == node(idx).clq_idx);
        clq_nbrs(find(clq_nbrs == idx)) = [];
        non_clq_nbrs = setdiff(node(idx).clq_connections,clq_nbrs);
        node(idx).ff_rec = non_clq_nbrs(:)';
    else
        node(idx).ff_rec = node(idx).clq_connections(:)';
    end
    ff_update = [ff_update idx];

end
end


