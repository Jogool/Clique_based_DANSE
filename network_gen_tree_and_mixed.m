% This function is used to generate a tree-topology for the T-DANSE
% algorithm as well as a clique based topology

function [node,uo,uo_mix] =  ... 
    network_gen_tree_and_mixed(nb_clq,min_clq_size,max_clq_size,nb_non_clq_nodes,dim_DANSE,nb_noise_srcs)

%close all
length_of_signal = 10000;
room_size = 5;
r_sen = 0.1;    % sensor distance from center of node (equispaced) now 10 centimeters
clq_dis = .5;   % distance clique nodes are spaced around the center
nb_sensors = 3;
circ_pos = linspace(0,2*pi,nb_sensors+1); % position of sensors around node

node_counter = 1;
clq_nodes = [];
for ii = 1:nb_clq
    clq_center(ii,:) =  room_size*rand(2,1);
    if eq(length(min_clq_size:max_clq_size),1)
        nb_clq_nodes = min_clq_size;
    else
        nb_clq_nodes = randsample(min_clq_size:max_clq_size,1);
    end
    index = node_counter : node_counter + nb_clq_nodes - 1;
    clique(ii).nodes = index;
    for jj = 1:nb_clq_nodes
        node(node_counter).pos = clq_center(ii,:)' -clq_dis+2*clq_dis*rand(2,1);
        node(node_counter).connections = index(find(node_counter ~= index));
        A_orig(node_counter,index(find(node_counter ~= index))) = 1;
        A_orig(index(find(node_counter ~= index)),node_counter) = 1;
        node(node_counter).isClique = 1;
        node(node_counter).clq_idx = ii;
        clq_nodes = [clq_nodes node_counter];
        node_counter  = node_counter + 1;
    end
end
non_clq_nodes = [];
for ii = 1:nb_non_clq_nodes
    node(node_counter).pos = room_size*rand(2,1);
    non_clq_nodes = [non_clq_nodes node_counter];
    node(node_counter).isClique = 0;
    node(node_counter).clq_idx = 0;
    node_counter  = node_counter + 1;
end

nb_nodes = size(node,2);
cmap = hsv(nb_nodes);
theta = linspace(0,2*pi,100);
% figure
% %axis([0 5 0 5]);
% hold on
% for ii = 1:nb_nodes
%     plot(node(ii).pos(1),node(ii).pos(2),'o','Color',cmap(ii,:));
%     %plot(node(ii).pos(1)+r*cos(theta),node(ii).pos(2)+r*sin(theta),'--','Color',cmap(ii,:))
%     text(node(ii).pos(1)+0.1,node(ii).pos(2)+0.1,num2str(ii))
% %     for jj = 1:numel(node(ii).connections)
% %         line([node(ii).pos(1) node(node(ii).connections(jj)).pos(1)],[node(ii).pos(2) node(node(ii).connections(jj)).pos(2)])
% %     end
% end
for ii = 1:nb_nodes
    index = find(ii ~= 1:nb_nodes);
    for jj = 1:nb_nodes-1
        D(ii,index(jj)) = norm(node(ii).pos-node(index(jj)).pos);
    end
end

% plot only center of cliques and non-clique nodes
% figure
% hold on
idx = [1:nb_clq non_clq_nodes];
pos_idx = [clq_center' [node(non_clq_nodes).pos]]';
for ii = 1:numel(idx)
%     plot(pos_idx(ii,1),pos_idx(ii,2),'o');
%     text(pos_idx(ii,1)+0.1,pos_idx(ii,2)+0.1,num2str(ii))
    temp_idx = find(ii ~= 1:numel(idx));
    for jj = temp_idx
        D_clq(ii,jj) = norm(pos_idx(ii,:)-pos_idx(jj,:));
    end
end
DG = sparse(D_clq);
[ST,~] = graphminspantree(DG);
A_tree = full(ST)+full(ST)';
A_tree(find(A_tree)) = 1;
%graphViz4Matlab(A_tree)

A_clq = zeros(nb_nodes);
A_clq(clq_nodes,clq_nodes) = A_orig(clq_nodes,clq_nodes);

%find shortest distnace between neighboring cliques
for ii = 1:nb_clq
    temp_idx = find(A_tree(ii,:));
    for jj = temp_idx
        if le(jj,nb_clq)
            [~,I] = min(reshape(D(clique(ii).nodes,clique(jj).nodes),length(clique(ii).nodes)*length(clique(jj).nodes),1));
            [I,J] = ind2sub([length(clique(ii).nodes),length(clique(jj).nodes)],I);
            A_clq(clique(ii).nodes(I),clique(jj).nodes(J)) = 1;
            A_clq(clique(jj).nodes(J),clique(ii).nodes(I)) = 1;
        end
    end
end
%graphViz4Matlab(A_clq)

% figure
% axis([0 5 0 5]);
% hold on
% for ii = 1:nb_nodes
%     plot(node(ii).pos(1),node(ii).pos(2),'o','Color',cmap(ii,:));
%     text(node(ii).pos(1)+0.1,node(ii).pos(2)+0.1,num2str(ii))
%     index = find(A_clq(ii,:));
%     for jj = index
%         line([node(ii).pos(1) node(jj).pos(1)],[node(ii).pos(2) node(jj).pos(2)])
%     end
% end

for ii = non_clq_nodes
    temp_idx = find(idx == ii);
    connections = find(A_tree(temp_idx,:)); % this gets all the connections connected to this node
    % if node is connected to clique, then find clique node with shortest
    % distance, else just make connection to non-clique node
    for jj = 1:numel(connections)
        if le(connections(jj),nb_clq)
            [~,I] = min(D(clique(connections(jj)).nodes,idx(temp_idx)));
            A_clq(clique(connections(jj)).nodes(I),idx(temp_idx)) = 1;
            A_clq(idx(temp_idx),clique(connections(jj)).nodes(I)) = 1;
        else
            A_clq(idx(temp_idx),idx(connections(jj))) = 1;
            A_clq(idx(connections(jj)),idx(temp_idx)) = 1;
        end
    end
end

% figure
% axis([0 5 0 5]);
% hold on
% for ii = 1:nb_nodes
%     plot(node(ii).pos(1),node(ii).pos(2),'o','Color',cmap(ii,:));
%     text(node(ii).pos(1)+0.1,node(ii).pos(2)+0.1,num2str(ii))
%     index = find(A_clq(ii,:));
%     for jj = index
%         line([node(ii).pos(1) node(jj).pos(1)],[node(ii).pos(2) node(jj).pos(2)])
%     end
% end
%graphViz4Matlab(A_clq)

%non_clq_nodes = find(cellfun(@isempty,{node.isClique}));


Deg_mat = diag(sum(A_clq,2));  % diagnomal matrix is equal to connections
L = Deg_mat-A_clq;
lambda = eig(L);
[~,I] = sort(lambda);
alg_conn_clq = lambda(I(2)); % algebratic connectivity

% rework tree so it now contains all nodes instead of just centers of
% cliques
DG = sparse(D);
[ST,~] = graphminspantree(DG);
A_tree = full(ST)+full(ST)';
A_tree(find(A_tree)) = 1;
%graphViz4Matlab(A_tree)

D_tree = diag(sum(A_tree,2));  % diagnomal matrix is equal to connections
L_tree = D_tree-A_tree;
lambda_tree = eig(L_tree);
[~,I] = sort(lambda_tree);
alg_conn_tree = lambda_tree(I(2));


% find the root node
[u1,v1]=eig(A_tree);
indicies = find(~A_tree);
[~,lambdaind]=sort(abs(diag(v1)));
EigvCentrality=[[1:size(A_tree,1)]' abs(u1(:,lambdaind(end)))];
[~,root_node] = max(EigvCentrality(:,2));   % pick root node based on centrality
node(root_node).isroot = 1;

% now generate signals
for ii = 1:dim_DANSE
    source(ii).nb = ii;
    source(ii).pos  = 5*rand(2,1);
    source(ii).signal = -0.5+rand(length_of_signal,1);
    
end
for ii = 1:nb_noise_srcs
    noise(ii).nb = ii;
    noise(ii).pos  = 5*rand(2,1);
    noise(ii).signal = -0.5+rand(length_of_signal,1);
end


white_noise_var = mean(var(cat(2,source.signal)))/2;


for ii = 1:nb_nodes
    node(ii).nb = ii;
    node(ii).sensors = nb_sensors;
    node(ii).ss_clean = zeros(length_of_signal,nb_sensors);
    node(ii).ss_noise = zeros(length_of_signal,nb_sensors);
    node(ii).sensors = nb_sensors;
    
    for jj = 1:nb_sensors
        node(ii).sensor(jj).pos(1) = node(ii).pos(1)+r_sen*cos(circ_pos(jj));
        node(ii).sensor(jj).pos(2) = node(ii).pos(2)+r_sen*sin(circ_pos(jj));
        for kk = 1: size(source,2)
            d = norm(source(kk).pos - node(ii).sensor(jj).pos');
            node(ii).steering(jj,kk) = d;
            node(ii).ss_clean(:,jj) = node(ii).ss_clean(:,jj)+d*source(kk).signal;
        end
        for kk = 1:size(noise,2)
            d = norm(noise(kk).pos - node(ii).sensor(jj).pos');
            node(ii).ss_noise(:,jj) = node(ii).ss_noise(:,jj)+d*noise(kk).signal+sqrt(white_noise_var)*randn(length_of_signal,1);
        end
    end
    
end
node = rmfield(node,'connections');
% now we have cliques connected a tree, now need to set up the tree using
% the same links

A_tree = A_clq;
for ii = 1:nb_clq
    clq_idx = find([node.clq_idx] == ii);
    A_tree(clq_idx,clq_idx) = 0;
    A_temp = A_clq;
    root_node_conn = find(A_temp(clq_idx,root_node));
    adj_power = 2;
    while isempty(root_node_conn)
        A_temp = A_clq^adj_power;
        root_node_conn = find(A_temp(clq_idx,root_node));
        adj_power = adj_power + 1;
    end
    if gt(numel(root_node_conn),1)
        A_tree(clq_idx(find(clq_idx ~= root_node)),clq_idx(find(clq_idx ~= root_node))) = 0;
        A_tree(clq_idx(find(clq_idx ~= root_node)),root_node) = 1;
        A_tree(root_node,clq_idx(find(clq_idx ~= root_node))) = 1;
       % A_tree(clq_idx(root_node_conn),clq_idx(root_node_conn)) = 0;
    else
        A_tree(clq_idx(find(clq_idx ~= clq_idx(root_node_conn))),clq_idx(root_node_conn)) = 1;
        A_tree(clq_idx(root_node_conn),clq_idx(find(clq_idx ~= clq_idx(root_node_conn)))) = 1;
    end
end
for ii = 1:nb_nodes
        node(ii).loc_filt_coeff = -1+2*rand(node(ii).sensors,dim_DANSE);
end

for ii = 1:nb_nodes
    node(ii).tree_connections = find(A_tree(:,ii));
    node(ii).clq_connections = find(A_clq(:,ii));
end

% uo = update order starting from root node
uo = root_node;
uo = path_find(uo,A_tree,root_node);
disp('');

uo_mix = root_node;
[node.df_update] = deal(0);
node(root_node).df_update = 1;
[uo_mix,~] = path_find_mixed(uo_mix,root_node,node);
disp('');
end

%% This is a recursive function to generate the update order of the tree
%% not needed for seqeuntial updating
function [uo] = path_find(uo,A,v)
%index = find(A(uo(end),:));
index = setdiff(find(A(v(end),:)),uo);
if isempty(index)
    uo = [uo uo(end-1)];
else
    for ii = 1:length(index)
        uo = path_find([uo index(ii)],A,index(ii));
    end
    disp('');
    if eq(v,uo(1))
        uo(end) = [];
    else
        uo = [uo uo(find(uo == v,1)-1)];
    end
end
end

function [uo,node] = path_find_mixed(uo,v,node)

if node(v).isClique
    clq_nodes = find([node.clq_idx] == node(v).clq_idx);
    nb_clq_nodes = numel(clq_nodes);
    index = clq_nodes(find(clq_nodes ~= v));
    index(find([node(index).df_update])) = [];
    
    only_clq = index(find(cellfun(@numel,{node(index).clq_connections}) == nb_clq_nodes - 1));
    [node(only_clq).df_update] = deal(1);
    uo = [uo only_clq];
    index(find(cellfun(@numel,{node(index).clq_connections}) == nb_clq_nodes - 1)) = [];
    
    for ii = 1:numel(index)
        if and(~eq(node(index(ii)).df_update,1),eq(numel(node(index(ii)).clq_connections),nb_clq_nodes - 1))
            uo = [uo index(ii)];
            node(index(ii)).df_update = 1;
        end
        if and(~eq(numel(node(index(ii)).clq_connections),nb_clq_nodes - 1),~eq(node(index(ii)).df_update,1))
            uo = [uo index(ii)];
            node(index(ii)).df_update = 1;
            [node(index).df_update] = deal(1);
            [uo,node] = path_find_mixed(uo,index(ii),node);
            [node(index).df_update] = deal(0);
            node(index(ii)).df_update = 1;
        end
    end
    if ~eq(v,uo(end))
        uo = [uo v];
        
    end
    non_clq_nodes = setdiff(node(v).clq_connections,clq_nodes);
    non_clq_nodes(find([node(non_clq_nodes).df_update])) = [];
    if ~isempty(non_clq_nodes)
        for ii = non_clq_nodes'
            uo = [uo ii];
            node(ii).df_update = 1;
            [uo,node] = path_find_mixed(uo,ii,node);
            if ~eq(v,uo(end))
                uo = [uo v];
            end
        end
    end
else
    index = node(v).clq_connections';
    index(find([node(index).df_update])) = [];
    for ii = index
        uo = [uo ii];
        node(ii).df_update = 1;
        [uo,node] = path_find_mixed(uo,ii,node);
    end
end
if and(node(v).isClique,isempty(index))
    disp('')
    [~,a,~] = intersect(uo,clq_nodes);
    a = sort(a);
    if eq(uo(a(1)),v)
        if eq(uo(1),uo(end))
            uo(end) = [];
        else
            uo = [uo uo(find(uo == v,1)-1)];
        end
    else
        uo = [uo uo(a(1))];
    end
    if eq(uo(end),v)
        uo(end) = [];
    end
else
    if eq(uo(1),uo(end))
        uo(end) = [];
    else
        uo = [uo uo(find(uo == v,1)-1)];
    end
end
disp('');

end

