% spantree   create a connected graph with N nodes
%
% SYNOPSIS:
%   edges = spantree(N)
%
% INPUT
%   N
%       The number of nodes
%
% OUTPUT
%   edges
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

function edges = spantree(N)

    nodeSet = 1:N;
    initIDX = randperm(N,1);
    initNode = nodeSet(initIDX);
    nodeSet(initIDX) = [];
    edges(1,1) = initNode;
    
    for i=1:N-1
        
        curIDX = randperm(N-i,1);
        curNode = nodeSet(curIDX);
        edges(2,i) = curNode;
        edges(1,i+1) = curNode;
        nodeSet(curIDX) = [];
        
    end       
    
    edges(:,N) = [];
    
end