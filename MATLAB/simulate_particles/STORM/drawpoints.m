function varargout = drawpoints(varargin)

% Check if input is a cell and reformat
if ((nargin == 2) && iscell(varargin{2}))
    if size(varargin{2},1) == 1
        varargin = cat(2,varargin{1},varargin{2});
    else
        varargin = cat(1,varargin{1},varargin{2});
    end
end 

% Read out variables
if numel(varargin) == 7
    particle_model = varargin{1};
    tech = varargin{2};
    timeframes = varargin{3};
    k_on = varargin{4};
    k_off = varargin{5};
    k_bleach = varargin{6};
    subframes = varargin{7};
else
    error('Wrong number of input arguments.')
end
    
% Initialize variables
N_mols = size(particle_model,1);            % Number of molecules/binding sites in the particle

% Simulate switching behavior
switch tech
    case 'palm'
        t = -1/k_on*log(1-rand(N_mols,1));                                                  % Localization times                                     
        locs_tmp = particle_model(t<=timeframes,:);                                         % Localized emitters in original reference frame
        t = t(t<=timeframes);
    case 'storm'
        % Times each molecule is localized
        M = min(binornd(timeframes,(1-exp(-k_on)),N_mols,1),1+geornd((1-exp(-k_bleach)),N_mols,1));   
        
        % Insert each nn-th particle position M(nn) times into locs_tmp
        if size(particle_model,2) == 2
            locs_tmp = zeros(sum(M),2);
        else
            locs_tmp = zeros(sum(M),3);
        end
        t = zeros(sum(M),1);
        st = 1;        
        for nn = 1:N_mols
            locs_tmp(st:st+M(nn)-1,:) = repmat(particle_model(nn,:),[M(nn) 1]);
            t(st:st+M(nn)-1) = randperm(timeframes,M(nn));%does not represent storm kinetics! looks like random over whole acquistion, flat histogram of #locs
            st = st + M(nn);
        end
        siteID = [];
        for i=1:N_mols
           siteID = [siteID;zeros(M(i),1)+i]; 
        end        
    case 'gsdim'
        locs_tmp = [];
        t = [];
        if subframes == 1
            active_mols = logical(binornd(1,k_on/(k_on+k_off),N_mols,1));
            bleached_mols = false(N_mols,1);
            
            % Simulate 
            for tt = 1:timeframes
                % Randomly activate initially inactive molecules
                active_mols(~active_mols) = logical(binornd(1,1-exp(-k_on),sum(~active_mols),1));
                
                % Add molecules that are active for some time
                locs_tmp = cat(1,locs_tmp,particle_model(active_mols,:));
                t = cat(1,t,tt*ones(sum(active_mols),1));

                % Randomly bleach molecules
                bleached_mols(active_mols) = logical(binornd(1,1-exp(-k_bleach),sum(active_mols),1));
                % Randomly leave initially active molecules active if they
                % did not bleach in between
                active_mols(active_mols) = (~bleached_mols(active_mols))&logical(binornd(1,exp(-k_off),sum(active_mols),1));
            end            
        else
            active_mols = logical(binornd(1,k_on/(k_on+k_off),N_mols,1));
            bleached_mols = false(N_mols,1);
            frac_active = [];
            % Simulate 
            for tt = 1:timeframes
                sub_active = zeros(N_mols,1);
                for ss = 1:subframes
                    % Randomly activate initially inactive molecules
                    active_mols(~active_mols) = logical(binornd(1,1-exp(-k_on/subframes),sum(~active_mols),1));
                    
                    % Count molecules that are active for some time
                    sub_active = sub_active+active_mols;

                    % Randomly bleach molecules
                    bleached_mols(active_mols) = logical(binornd(1,1-exp(-k_bleach/subframes),sum(active_mols),1));
                    % Randomly leave initially active molecules active if they
                    % did not bleach in between
                    active_mols(active_mols) = (~bleached_mols(active_mols))&logical(binornd(1,exp(-k_off/subframes),sum(active_mols),1));
                end
               
                % Add molecules that are active for some time
                locs_tmp = cat(1,locs_tmp,particle_model(sub_active>0,:));
                t = cat(1,t,tt*ones(sum(sub_active>0),1));
                frac_active = cat(1,frac_active,sub_active(sub_active>0)/subframes);
            end 
            varargout{3} = frac_active;
        end
        
    otherwise
        locs_tmp = particle_model;
end

% Write outputs
varargout{1} = locs_tmp;
varargout{2} = t;

if exist('frac_active','var')
    varargout{3} = frac_active;
else
    varargout{3} = ones(size(t));
end
varargout{4} = siteID;