path_matlab = genpath('Optimization');
addpath(path_matlab)
cvx_solver mosek;
cvx_solver

load('../data/simulation_final/tud_nph5000_dol100_tr10_nmN256sym.mat')
% comment
% Particles = Particles(1:128);

mainroot = '../';
allFiles = dir( mainroot );
allFiles(1:2) = [];
nMachine = numel(allFiles);
    
RM = zeros(4,4,5);
I = zeros(2,5);
iter = 1;

    for i=1:1%nMachine
       
        allRows = dir([mainroot 'test_tuddol100sim/']);
        allNames = {allRows(~[allRows.isdir]).name};
        allNames = natsortfiles(allNames);
        nRows = numel(allNames);
        
        for j=1:nRows
            
            load([mainroot 'test_tuddol100sim/' allNames{1,j}]);
            N = numel(result);
                     
            for k=1:N
                
                I(:,iter) = result(k).id;

                estAngleR = -result(k).parameter(3);
                curestAngle = wrapToPi(estAngleR);
                w_est = curestAngle * [0;0;1];
                RM(1:3,1:3,iter) = expm([0,-w_est(3),w_est(2); w_est(3),0,-w_est(1); -w_est(2),w_est(1),0]);
                RM(:,4,iter) = [result(k).parameter(1:2)'; 0; 1];                
               
                iter = iter + 1;
            end
            
        end
        
    end

% uncomment for experimental dataset
% for i=1:numel(particles)
%     Particles{1,i}.coords(:,1:2) = particles{1,i}.points;
% end
% for i=276:549
%     Particles{1,i-275}.coords(:,1:2) = particles{1,i}.points;
% end


nParticles = 256;%numel(Particles);

totalPairs = 0.5 * nParticles * (nParticles-1);

edges = motion_graph(nParticles);
for i=1:nParticles-1
    idxPair(i) = find(I(1,:) == edges(1,i) & I(2,:) == edges(2,i));
end
nPairs = totalPairs - numel(idxPair); % change to average only a subset
totalPairsIDX = 1:totalPairs;
totalPairsIDX(idxPair) = [];
idxRest = randperm(numel(totalPairsIDX), nPairs);
idxPair = [idxPair totalPairsIDX(idxRest)];
[M] = MeanSE3Graph(RM(:,:,idxPair),I(:,idxPair));

superParticle0 = [];
for i = 1:nParticles%numel(Particles)
        
     par = Particles{1,i}.coords(:,1:2);
     tmpParticle = (par -repmat(M(1:2,4,i)',size(par,1),1))* M(1:2,1:2,i)';
     superParticle0 = [superParticle0; tmpParticle];

end
%%
relTr_init = zeros(2,5);
% initial relative all2all matrix
for i=1:totalPairs
    relAngleR_init(i) = atan2(RM(2,1,i),RM(1,1,i));
    relTr_init(:,i) = RM(1:2,4,i);
end

% relative all2all matrix after averaging
kk = 1;
relTr = zeros(2,5);
for i=1:nParticles-1
    for j=i+1:nParticles
        
        relAngleR(kk) = atan2(M(2,1,i),M(1,1,i)) - atan2(M(2,1,j),M(1,1,j));
        relTr(:,kk) = M(1:2,4,i)-M(1:2,4,j);
        kk = kk+1;
        
    end
end
%%
error = wrapTo180(rad2deg(relAngleR + relAngleR_init));
error_tr = sqrt(sum((relTr - relTr_init).^2,1));
error_idx = find(abs(error) > 5);
% error_idx = find(abs(error) > 5 | error_tr > (0.025));
%%
% RM_new = M2RM(M);
RM_new = RM;
I_new = I;

RM_new(:,:,error_idx) = [];
I_new(:,error_idx) = [];

[M_new] = MeanSE3Graph(RM_new,I_new);

superParticle1 = [];
for i = 1:nParticles%numel(Particles)
        
     par = Particles{1,i}.coords(:,1:2);
     tmpParticle = (par -repmat(M_new(1:2,4,i)',size(par,1),1))* M_new(1:2,1:2,i)';
     superParticle1 = [superParticle1; tmpParticle];

end

%%
%{
% relative all2all matrix after averaging
kk = 1;
for i=1:nParticles-1
    for j=i+1:nParticles
        
        relAngleR2(kk) = atan2(M_new(2,1,i),M_new(1,1,i)) - atan2(M_new(2,1,j),M_new(1,1,j));
        kk = kk+1;
        
    end
end

error2 = wrapTo180(rad2deg(relAngleR2 - relAngleR));
% error_idx2 = find(abs(error2) > 5);
error_idx2 = find(error_tr > (0.02));

RM_new2 = RM;
I_new2 = I;

RM_new2(:,:,error_idx2) = [];
I_new2(:,error_idx2) = [];

[M_new2] = MeanSE3Graph(RM_new2,I_new2);

superParticle4 = [];
for i = 1:numel(Particles)
        
     par = Particles{1,i}.coords(:,1:2);
     tmpParticle = (par -repmat(M_new2(1:2,4,i)',size(par,1),1))* M_new2(1:2,1:2,i)';
     superParticle4 = [superParticle4; tmpParticle];

end
%}
%%
save(['simulation_final/dol100test/motion_' num2str(nPairs)], 'M');
save(['simulation_final/dol100test/particle_' num2str(nPairs)], 'superParticle0');

save(['simulation_final/dol100test/motion_outlier_removal_' num2str(nPairs)], 'M_new');
save(['simulation_final/dol100test/particle_outlier_removal_' num2str(nPairs)], 'superParticle1');

% save(['result_Motion_new_new' num2str(nPairs)], 'M_new2');
% save(['result_particle_new_new' num2str(nPairs)], 'superParticle3');
