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
DANSE_param.nb_nodes = 9;       

% number of sensors per node (assumed same across all nodes compression 
%ratio of (DANSE_param.sensors+1)/DANSE_param.desired_sources)
DANSE_param.sensors = DANSE_param.desired_sources + 1; 

% number of correlated noise sources
DANSE_param.noise_sources = 4;    

% number and size of cliques 
DANSE_param.nb_clq = 2;     % number of cliques
DANSE_param.clq_size = 3;   % number of nodes in clique
if lt(DANSE_param.nb_nodes,DANSE_param.nb_clq*DANSE_param.clq_size)
    disp(['Warning not enough nodes for clique generation']);
end

plot_on = 1;        % 1(0) - show (do not show) network

% max iterations before stopping DANSE algorithms 
%note algoritm may not have converged when max iter has been reached
max_iter = 500;  

 % threshold for when to stop algorithms, i.e., when convergence is met
thresh = 1e-4;     

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

% remove redudant updating from clique path : removes all clique only 
% nodes (only let them update once)

updateorder_clq = updateorder;
idx = find([node.isclq]);
for ii = 1:numel(idx) 
    if isempty(setdiff(node(idx(ii)).clq_conn,node(idx(ii)).clq_nbrs))
        [~,I] = find(updateorder_clq == idx(ii));
        updateorder_clq(I(2:end)) = [];
    else
        [~,I] = find(updateorder_clq == idx(ii));
        idx_rem = [];
        for jj = 1:numel(I)
            % chech to make sure that node is not at the end or the
            % beginning of the update path
            if and(~eq(I(jj),1),lt(I,numel(updateorder_clq)))
                start_node = updateorder_clq(I(jj)-1);
                end_node = updateorder_clq(I(jj)+1);
                [~,~,ib] = intersect(node(idx(ii)).clq_nbrs,[start_node end_node]);
                if eq(numel(ib),2)
                    idx_rem = [idx_rem I(jj)];
                end
            end
        end
        updateorder_clq(idx_rem) = [];
    end
end

% find centralized solution
disp('Centralized')
[node] = centralized(node);
%load node
save node node updateorder_clq updateorder
savefig('node.fig')
% store original coefficients, this is loaded before every instance of a
% DANSE algorithm, so that the local filters all start at the same value
org_node = node;

%% DANSE - round robin updating
disp('DANSE round robin updating')
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

%% DANSE - tree updating
node = org_node;
fprintf('\n')
disp('DANSE tree updating')
reverseStr = '';

node_update = 1;
cost_sum_DANSE_tree_updating = [];
ii = 1;
while 1
    [node] = DANSE(node,node_update);
    cost_sum_DANSE_tree_updating = [cost_sum_DANSE_tree_updating ...
        sum(cat(1,node.cost))];
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
%% T-DANSE
node = org_node;
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

%% DANSE mixed topology - clique based updating order
node = org_node;
fprintf('\n')
reverseStr = '';
disp('MTDANSE MUO')
node_update = updateorder(1);
cost_sum_MTDANSE_MUO = [];
ii = 1;
while 1
    if or(eq(ii,32),or(eq(ii,46),eq(ii,60)))
        disp('')
    end
    [node] = MTDANSE(node,node_update);
    cost_sum_MTDANSE_MUO = [cost_sum_MTDANSE_MUO sum(cat(1,node.cost))];
    tot_diff = norm(cat(1,node.cost_cent) - ...
        cellfun(@(x) x(end), {node.cost})');
    
    if or(lt(tot_diff,thresh),ge(ii,max_iter));
        break
    else
        ii = ii + 1;
    end
    node_update=updateorder_clq(rem(ii,numel(updateorder_clq))+1);
    msg = sprintf('Iteration : %d', ii);
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
end
fprintf('\n')
%% DANSE mixed topology - tree based updating order
node = org_node;
reverseStr = '';
disp('MTDANSE TUO')
node_update = updateorder(1);
cost_sum_MTDANSE_TUO = [];
ii = 1;
while 1
    [node] = MTDANSE(node,node_update);
    cost_sum_MTDANSE_TUO = [cost_sum_MTDANSE_TUO sum(cat(1,node.cost))];
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
fprintf('\n')

%% plot
if plot_on
    figure
    hold on
    loglog(cost_sum_DANSE)
    loglog(cost_sum_DANSE_tree_updating,'-xm')
    loglog(cost_sum_TDANSE,'-or')
    loglog(cost_sum_MTDANSE_MUO,'--dk')
    loglog(cost_sum_MTDANSE_TUO,'-*g')
    axis tight

h =  refline(0,sum([node.cost_cent]));
set(h,'LineStyle','--');

a = get(gca,'YLim');
set(gca,'YLim',[sum([node.cost_cent]) - sum([node.cost_cent])*.1 a(2)])
legend('DANSE - FC-RR', 'DANSE - FC-TUO', 'T-DANSE','DANSE - MT-MUO','DANSE - MT-TUO', 'Optimal');
set(gca, 'XScale', 'log', 'YScale', 'log')
xlabel('Iteration')
ylabel('Sum of LS cost for all nodes (dB)')
end

%total_conv(1,1:length(cost_sum_DANSE)) = cost_sum_DANSE;
%total_conv(2,1:length(cost_sum_TDANSE)) = cost_sum_TDANSE;
%total_conv(3,1:length(cost_sum_TIDANSE_fc)) = cost_sum_TIDANSE_fc;
%total_conv(4,1:length(cost_sum_TIDANSE_tree)) = cost_sum_TIDANSE_tree;
%total_conv(5,1) = sum([node.cost_cent]);



