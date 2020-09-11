% Cell_out Convert a structure with function parameters to a cell array
%
% An input structure with parameter definitions is provided that is structured
% according to the conventions of the inparams stucture for getparams.
% The other input structure is searched for fields with these parameters
% and put into a cell array according to the ordering in the definitions.

% SYNOPSIS:
%   function cell_out = convertparamstruct(defs,paramstruct_in)
%
% (C) Copyright 2012               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
% Robert Nieuwenhuizen, Jan 2013
%
function out = convertparamstruct(varargin)

% Check number of inputs
if nargin ~= 2
    error('Too few input arguments given.');
end

% Check if inputs are structures
if ~(isstruct(varargin{1}))
    error('No valid definitions provided.');
end

if ~isstruct(varargin{2})
    error('No parameters for conversion specified.');
end

% Initialize variables
defs = varargin{1};
paramstruct_in = varargin{2};
command_string = '[';
paramcell = cell(1,numel(defs));

% Add the variables from defs that are in paramstruct_in to paramcell in
% the order specified by defs
for ii = 1:numel(defs)
    command_string = [command_string,defs(ii).name,','];        % Write variable name to string for getparams

    if isfield(paramstruct_in,defs(ii).name)
        paramcell{ii} = getfield(paramstruct_in,defs(ii).name);         
    else
        paramcell{ii} = defs(ii).default;
    end        
end

command_string = [command_string(1:end-1),']'];
defs = struct('inparams',defs);

% Check if provided variables satisfy the criteria specified in defs
try
    eval([command_string,' = getparams(defs,paramcell{:});']);
catch
   if ~isempty(paramerror)
      error(paramerror)
   else
      error('Parameter parsing was unsuccessful.')
   end
end

% Output arguments as a cell array for easy passing to functions
out = eval(['{',command_string(2:end-1),'}']);