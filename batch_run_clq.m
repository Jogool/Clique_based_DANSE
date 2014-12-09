function [total_conv] = batch_run_clq()
% Syntax:  [cout] = batch_run()
%
% Inputs:   none - for options change hardcoded parameters
%                                                         
%
% Outputs: total_conv       - store the summed LS cost of all nodes
%                             for all algorithms

% Other m-files required: network_gen_tree_clq, centralized, folders
% (DANSE,TDANSE,DANSE_clq)
% MAT-files required: none
%
% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% Dec. 2014; Last revision: 09-Dec-2014

%% hardcoded parameters
% number of desired sources (also dimension of DANSE)
DANSE_param.desired_sources = 2;  

% number of nodes
DANSE_param.nb_nodes = 15;       

% number of sensors per node (assumed same across all nodes compression 
%ratio of (DANSE_param.sensors+1)/DANSE_param.desired_sources)
DANSE_param.sensors = DANSE_param.desired_sources + 1; 

% number of correlated noise sources
DANSE_param.noise_sources = 4;    

% number and size of cliques 
DANSE_param.nb_clq = 3;     % number of cliques
DANSE_param.clq_size = 4;   % number of nodes in clique
if lt(DANSE_param.nb_nodes,DANSE_param.nb_clq*DANSE_param.clq_size)
    disp(['Warning not enough nodes for clique generation']);
end

plot_on = 1;        % 1(0) - show (do not show) network

% max iterations before stopping DANSE algorithms 
%note algoritm may not have converged when max iter has been reached
max_iter = 1000;  

 % threshold for when to stop algorithms, i.e., when convergence is met
thresh = 1e-5;     

% output 
total_conv = zeros(5,max_iter); % see header for description

% generate random network and TDANSE updating order 
if plot_on
    close all
    [node,source,noise,wnv] = network_gen_clq(DANSE_param);
    [node,updateorder] = construct_tree_clq(node);
    plot_WSN_clq(node,source,noise)
else
    [node,~,~,~] = network_gen_clq(DANSE_param);
    [node,~,updateorder] = construct_tree_clq(node);
end

% remove redudant updating from clique path
% this removes all clique only nodes, only let them update once
% something is wrong with removing only clique connections 
updateorder_clq = updateorder;
idx = find([node.isclq]);
idx = setdiff(idx',updateorder_clq');
for ii = 1:numel(idx)
    if eq(numel(node(idx(ii)).clq_conn),numel(node(idx(ii)).clq_nbrs))
        [~,~,ib] = intersect(node(idx(ii)).clq_nbrs,updateorder_clq);
        ib = sort(ib);
        ib = ib(1);
        updateorder_clq = [updateorder_clq(1:ib) idx(ii) updateorder_clq(ib+1:end)];
        %if eq(ib(1)+1,ib(2))
        %    updateorder_clq = [] ;
        %else
            
        %end
    end
end

[node(updateorder_clq).clq];
updateorder_clq
%idxs = [node(updateorder_clq).clq];
idx = 1;
counter = 0;
idx_start = 0;
uo_temp = [];
while lt(idx,numel(updateorder_clq))
    if eq(idxs(idx_start+1),idxs(idx));
        counter = counter + 1;
    else
        if gt(counter,DANSE_param.clq_size)
            temp_idx_start = updateorder_clq(idx_start+1);
            temp_idx_end = updateorder_clq(idx-1);
            clq_idx = node(temp_idx_start).clq_nbrs;
            clq_idx(find(clq_idx == temp_idx_end)) = [];
            updateorder_temp = [temp_idx_start  clq_idx temp_idx_end];
            uo_temp = [uo_temp updateorder_temp];
            
        else
            if ~idx_start
                idx_start = 1;
            end
            uo_temp = [uo_temp updateorder_clq(idx_start:idx-1)];
        end
        counter = 0;
        idx_start = idx;
    end
    idx = idx + 1;
end
updateorder_clq
uo_temp

% find centralized solution
disp('Centralized')
[node] = centralized(node);

% store original coefficients, this is loaded before every instance of a
% DANSE algorithm, so that the local filters all start at the same value
org_node = node;

% DANSE
fprintf('\n')
disp('DANSE')
reverseStr = '';

node_update = 1;
cost_sum_DANSE = [];
ii = 1;
while 1
    [node] = DANSE(node,node_update);
    cost_sum_DANSE = [cost_sum_DANSE sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) - ...
        cellfun(@(x) x(end), {node.cost})');
    
    if or(lt(tot_diff,thresh),ge(ii,max_iter));
        break
    else
        ii = ii + 1;  
    end
    node_update=rem(node_update,DANSE_param.nb_nodes)+1;    
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

node = org_node;
% T-DANSE
fprintf('\n')
reverseStr = '';
disp('TDANSE')
node_update = updateorder(1);
cost_sum_TDANSE = [];
ii = 1;
while 1
    [node] = TDANSE(node,node_update);
    cost_sum_TDANSE = [cost_sum_TDANSE sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) - ...
        cellfun(@(x) x(end), {node.cost})');
    
    if or(lt(tot_diff,thresh),ge(ii,max_iter));
        break
    else
        ii = ii + 1;
    end
    node_update=updateorder(rem(ii,numel(updateorder))+1);
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end

node = org_node;
fprintf('\n')
reverseStr = '';
disp('Clique based')


if plot_on
    figure
    hold on
    loglog(cost_sum_DANSE)
    loglog(cost_sum_TDANSE,'-xm')
    loglog(cost_sum_TIDANSE_fc,'-or')
    loglog(cost_sum_TIDANSE_tree,'--dk')
    axis tight

h =  refline(0,sum([node.cost_cent]));
set(h,'LineStyle','--');

a = get(gca,'YLim');
set(gca,'YLim',[sum([node.cost_cent]) - sum([node.cost_cent])*.1 a(2)])
legend('DANSE', 'T-DANSE', 'TI-DANSE (FC)','TI-DANSE (T)', 'Optimal');
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Iteration')
ylabel('Sum of LS cost for all nodes (dB)')
end

total_conv(1,1:length(cost_sum_DANSE)) = cost_sum_DANSE;
total_conv(2,1:length(cost_sum_TDANSE)) = cost_sum_TDANSE;
total_conv(3,1:length(cost_sum_TIDANSE_fc)) = cost_sum_TIDANSE_fc;
total_conv(4,1:length(cost_sum_TIDANSE_tree)) = cost_sum_TIDANSE_tree;
total_conv(5,1) = sum([node.cost_cent]);
fprintf('\n')



