% Simlocs_particles Simulates a 3d acquisition of single particles with random orientations
%
% function [positions, angles] = simlocs_particles(particle_shape,N,acquisition_params,acquisition_3d)
%
function varargout = simlocs(varargin)

d = struct('menu','FRC resolution',...
           'display','Single particle simulation',...
           'inparams',struct('name',       {'particle_model',    'N',                   'switching_params',     'localization_params'},...
                             'description',{'Particle model',    'Number of particles', 'Switching parameters', 'Localization parameters'},...
                             'type',       {'array',             'array',               'cellarray',            'cellarray'},...
                             'dim_check',  {[-1 2],              0                      [],                     []},...
                             'range_check',{[],                  'N+',                  [],                     []},...
                             'required',   {1,                   1,                     0,                      0},...
                             'default',    {[],                  1,                     {[]},                   {[]}}...
                              ),...
           'outparams',struct('name',{'positions',              'angles'},...
                              'description',{'Localizations',   'True Euler angles'},...
                              'type',{'array',                  'array'}...
                              )...
           );       
       
if nargin == 1
   s = varargin{1};
   if ischar(s) && strcmp(s,'DIP_GetParamList')
      varargout{1} = struct('menu','none');
      return
   end
end

try
    switching_params = varargin{3};                         % See below for further parsing
    localization_params = varargin{4};
    [particle_model, N, ~, ~] = getparams(d,varargin{1:2});
catch
    if ~(exist('switching_params','var') && exist('localization_params','var'))
        error('Wrong number of input variables.');
    else
        if ~isempty(paramerror)
            error(paramerror)
        else
            error('Parameter parsing was unsuccessful.')
            
        end
    end
end

%% Read inputs

% Initialize variables
positions = [];

% Convert switching_params to a cell array with elements specified in
% switching_defs
defs_switching = struct('name',       {'tech',                    'timeframes',   'k_on',     'k_off',    'k_bleach',     'subframes'},...
                        'description',{'Switching mechanism',     'Time frames',  'On-rate',  'Off-rate', 'Bleach rate',  'Subframes'},...
                        'type',       {'string',                  'array',        'array',    'array',    'array',        'array'},...
                        'dim_check',  {0,                         0,              0,          {[],0},     0,              {[],0}},...
                        'range_check',{{'palm','storm'},          'R+',           'R+',       'R+',       'R+',           'N+'},...
                        'required',   {1,                         1,              1,          0,          0,              0},...
                        'default',    {'palm',                    1,              1,          eps,          eps,            1}...
                        );
switching_params = convertparamstruct(defs_switching,switching_params);

% Read out localization parameters
Ephotons = localization_params.photons;
bg = localization_params.bg;
PSFsigma_xy = localization_params.PSFsigma_xy;

%% Simulate localizations
% Draw points from particle model
[locs_tmp, t, ~, siteID] = drawpoints(particle_model,switching_params);   

%% Outputs
varargout{1} = locs_tmp;               % [x y]
varargout{2} = t;                      % [t]
varargout{3} = siteID;
