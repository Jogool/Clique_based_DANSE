function [node] = centralized(node)
%centralized - Finds the centralized solution 
%
% The function first generates a set of nodes randomly placed in a sensing
% environment.  The communication range of each node is increased
% until the network is connected.  Using this ad-hoc connectivity, a
% minimum spanning tree is then found.
%
% Syntax:  [node] = centralized(node)
%
% Inputs:
%   node            -   node structure from network_gen_tree.m
%
% Outputs:
%   node            -   node structure that contains the centralized 
%                       cost at each node
%
% Example: 
%    [node] = centralized(node)
%
% Other m-files required: none
% Subfunctions: node
% MAT-files required: none
%

% Author: Joseph Szurley
% email: joseph.szurley@esat.kuleuven.be
% March 2014; Last revision: 05-March-2004

%   assume nodes have access to all sensor signals

% gather all signals
x = [node.ss_clean];
n = [node.ss_noise];

% generate correlation matricies
Rnn = n'*n;
Rxx = x'*x;
y = x + n;
index = 1;
dim_DANSE = node(1).dimDANSE;
for ii = 1:size(node,2)
    %   calculate centralized filter
    w_temp = (Rnn+Rxx) \ Rxx(:,index:index+dim_DANSE-1);

    %   calculate centralized cost
    node(ii).cost_cent = norm(x(:,index:index+dim_DANSE-1)' - w_temp'*y')^2;
    %   store centralized local filter coefficents
    node(ii).local_filter_cent = w_temp(index:index+node(ii).sensors-1,:);
    node(ii).cent_filt = w_temp;
    index = index+node(ii).sensors;
end