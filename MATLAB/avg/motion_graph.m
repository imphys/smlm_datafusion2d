% motion_graph   create a connected graph with N nodes
%
% SYNOPSIS:
%   idx = motion_graph(N)
%
% INPUT
%   N
%       The number of nodes
%
% OUTPUT
%   idx
%       The edges of the graph, a 2x(N-1) matrix with the node indices as
%       its elements
%
% NOTES
%   
%
% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
%
% Hamidreza Heydarian, 2017

function idx = motion_graph(N)

    % initial minimum spanning tree
    for i=1:N-1
        spanning_tree(1,i) = i;
        spanning_tree(2,i) = i+1;
    end
    
    idx = spanning_tree;
    
end