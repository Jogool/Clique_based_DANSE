function [uo] = path_find_clq(uo,A,v,node)
% function to find updating order of TDANSE algorithm
%
% This is a recursive function to generate the updaing order of the TDANSE
% algorithm
%
% Syntax:  [node] = path_find(uo,A,v)
%
% Inputs:
%   uo          -   current updating order
%   A           -   adjacency matrix of tree
%   v           -   current node
%
% Outputs:

%   uo          -   updating order including next node
%                   
% Example:
%    [uo] = path_find(uo,A,v)
%

% Author: Joseph Szurley
% Work address
% email: joseph.szurley@esat.kuleuven.be
% December 2014; Last revision: 04-Dec-2014
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