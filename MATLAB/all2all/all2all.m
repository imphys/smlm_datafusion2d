% all2all   performs all 2 all registration for a given set of particles
%
% SYNOPSIS:
%   all2all(Particles, outdir)
%
% INPUT
%   Particles
%       Cell arrays of particles with localization in the point field and 
%       uncertainties in the sigma field.
%   outdir
%       Output directory where rows of all2all matrix are stored
%   scale
%       scale parameter for gmm registration
%
% OUTPUT
%   The function doesn't return any output but the results are stored in
%   outdir. Each element of the all2all matrix includes registration
%   parameters (result.parameter), cost function value (result.val) and the
%   indicator pair (result.id) which stores (row, column) indices of all2all 
%   matrix
%
% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
%
% Hamidreza Heydarian, 2017

function all2all(Particles, outdir, scale)

    % setup pyramid, determine the pyramid height or the number of layers
    N = numel(Particles);
    
    disp('all2all registration started !');
    disp(['There are ' num2str(N-1) ' rows !']);
    for idxM = 1:N-1
        tic
        disp(['row ' num2str(idxM) ' started!']);
        parfor idxS = idxM+1:N

            M = Particles{idxM};
            S = Particles{idxS};
            
            % perform pairwise registration for each element of all2all
            % matrix
            [param, ~, ~, ~, val] = pairFitting(M, S, scale);

            % registration parameters, cost function value and idicators
            % are stored in the result structure
            result(idxS-idxM).parameter = param;
            result(idxS-idxM).val = val;
            result(idxS-idxM).id = [idxM;idxS];

        end

%       save each row of the all2all matrix in result_xx where xx is the row number 
        if ~exist(outdir,'dir')
            mkdir(outdir);
        end        
        save([outdir '/result_' num2str(idxM)], 'result');
        clearvars result;
        a = toc;
        disp(['row ' num2str(idxM) ' done in ' num2str(a) ' seconds']);

    end

    disp('all2all registration done !');
    fprintf('\n\n');
    
end







