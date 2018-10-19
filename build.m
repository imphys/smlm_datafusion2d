% this script compiles the indicated source files below as well as related
% headers to produce a .mex binary which can be executed by Matlab.
% the resulting functionality is somewhat analogous to the CMakeLists.txt

clear all; %#ok<CLALL>

% "fullfile" pieces together strings into a full path
% "pwd" (obviously) returns a string of the working directory path
% "fileparts(pwd)" returns a string of the path one directory up

% include path
include_g = ['-I' ];

% specify source files
src_g_cpu = fullfile(pwd, 'gausstransform/mex_gausstransform_cpu.cpp');
src_g_gpu = fullfile(pwd, 'gausstransform/mex_gausstransform.cpp');
src_e_cpu = fullfile(pwd, 'expdist/mex_expdist_cpu.cpp');
src_e_gpu = fullfile(pwd, 'expdist/mex_expdist.cpp');

current_folder = pwd

% build
mex('-v', include_g, src_g_cpu);
mex('-v', include_e, src_e_cpu);

prompt = 'Also build mex files using GPU code (yes/no)? ';
gpu = input(prompt)

if gpu == 'yes'
    mex('-v', include_g, src_g_gpu, '-L', current_folder, '-l', 'gausstransform');
    mex('-v', include_e, src_e_gpu, '-L', current_folder, '-l', 'expdist');
end
