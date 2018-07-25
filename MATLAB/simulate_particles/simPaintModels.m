% simPaint   simulate SMLM particles from the model in simPaintModels
%
% SYNOPSIS:
%   [ Model ] = simPaintModels( model_name, raster_spacing, CCD_pixelsize )
%
% INPUT
%     model_name        the model name for SMLM log
%     raster_spacing    [nm]
%     CCD_pixelsize     [nm]
%
% OUTPUT
%   Model 
%       coordinates of the binding sites
%
% NOTES
% 
% (C) Copyright 2017               Quantitative Imaging Group
%     All rights reserved          Faculty of Applied Physics
%                                  Delft University of Technology
%                                  Lorentzweg 1
%                                  2628 CJ Delft
%                                  The Netherlands
%
% Author: Robert Nieuwenhuizen, Emile Boerder
% Edit:   Hamidreza Heydarian, 2017

function [ Model ] = simPaintModels( model_name, raster_spacing, CCD_pixelsize )

switch model_name
    
    case 'grid'
        Model = model_GRID(raster_spacing, CCD_pixelsize);
    
    case 'tud'
        Model = model_TUD(raster_spacing, CCD_pixelsize);
        
    case 'tud_ralf'
        Model = model_TUD_Ralf(raster_spacing, CCD_pixelsize);
        
    case 'x'
        Model = model_X(raster_spacing, CCD_pixelsize);

    case 'r'
        Model = model_R(raster_spacing, CCD_pixelsize);
        
    case 'hello'
        Model = model_HELLO(raster_spacing, CCD_pixelsize);

    case 'hello world'
        Model = model_HELLOWORLD(raster_spacing, CCD_pixelsize);

    case 'mpi'
        Model = model_MPI(raster_spacing, CCD_pixelsize);
        
    case 'phantum'
        Model = model_phantum(raster_spacing, CCD_pixelsize);  
        
    case 'phantum_wide'
        Model = model_phantum_wide(raster_spacing, CCD_pixelsize);  
end

%Model = Model + randn(size(Model)).*(raster_spacing/10)./CCD_pixelsize;

end

function [ Model ] = model_GRID( raster_spacing, CCD_pixelsize )

        s = 5;
        
        Model = zeros(s^2, 2);
        for i = 1:s
            Model(s*(i-1) + (1:s), 1) = i;
            Model(s*(i-1) + (1:s), 2) = 1:s;
        end

        Model = Model*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end

function [ Model ] = model_TUD( raster_spacing, CCD_pixelsize )

        %TUD
        Model = [

        %FLAME
        [-0.5      1];
        [-1.5      1];
        [-2.5      1];
        [-2      2];
        [-3      2];
        [-0.5      3];
        [-1.5      3];
        [-2.5      3];
        [-3.5      3];
        [-4      4];

        %T
        [1      0];
        [1.5    1];
        [1      2];
        [2      2];
        [3      2];
        [4      2];
        [5      2];
        [1.5    3];
        [1      4];

        %U
        [1      6];
        [2      6];
        [3      6];
        [4      6];
        [5      6];
        [5.5    7];
        [5.5    8];
        [1.5      9];
        [2.5      9];
        [3.5      9];
        [4.5      9];

        %D
        [1.5      11];
        [2.5      11];
        [3.5      11];
        [4.5      11];
        [5.5      11];
        [1      12];
        [5.5      12];
        [1.5      13];
        [5.5      13];
        [2.5      14];
        [3.5      14];
        [4.5      14];

        ]*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end

function [ Model ] = model_TUD_Ralf( raster_spacing, CCD_pixelsize )

        %TUD
        Model = [

        %FLAME
        [10      30];
        [10      35];
        [15      37.5];
        [20      30];
        [20      35];
        [20      40];
        [25      42.5];

        %T
        [0      20];
        [5      22.5];
        [10     5];
        [10     10];
        [10     15];
        [10     20];
        [15     22.5];
        [20     20];

        %U
        [25      7.5];
        [25     12.5];
        [25     17.5];
        [25     22.5];
        [30      5];
        [35      2.5];
        [40      5];
        [45      7.5];
        [45     12.5];
        [45     17.5];
        [45     22.5];

        %D
        [55      7.5];
        [55     12.5];
        [55     17.5];
        [55     22.5];
        [60      5];
        [65      7.5];
        [70      10];
        [70      15];
        [70      20];
        [60      25];
        [65      22.5];

        ]*0.2*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end

function [ Model ] = model_X( raster_spacing, CCD_pixelsize )

        %X
        Model = [
        
        [-8     -4];
        [-8     8];
        [-7     -3.5];
        [-7     7.5];
        [-6     -3];
        [-6     7];
        [-5     -2.5];
        [-5     6.5];
        [-4     -2];
        [-4     6];
        [-3     -1.5];
        [-3     5.5];
        [-2     -1];
        [-2     5];
        [-1     -0.5];
        [-1     4.5];
        [0      0];
        [0      4];
        [1      0.5];
        [1      3.5];
        [2      1];
        [2      3];
        [3      1.5];
        [3      2.5];
        [4      2];
        [16     -4];
        [16     8];
        [15     -3.5];
        [15     7.5];
        [14     -3];
        [14     7];
        [13     -2.5];
        [13     6.5];
        [12     -2];
        [12     6];
        [11     -1.5];
        [11     5.5];
        [10     -1];
        [10     5];
        [9     -0.5];
        [9     4.5];
        [8      0];
        [8      4];
        [7      0.5];
        [7      3.5];
        [6      1];
        [6      3];
        [5      1.5];
        [5      2.5];

        ]*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end

function [ Model ] = model_HELLO( raster_spacing, CCD_pixelsize )

        %HELLO
        Model = [

        %H
        [1      0];
        [1      1];
        [1      2];
        [1      3];
        [1      4];
        [1      5];
        [1      6];
        [2      0.5];
        [2      1.5];
        [2      2.5];
        [2      3.5];
        [2      4.5];
        [2      5.5];
        [3      2];
        [3      3];
        [4      2.5];
        [4      3.5];
        [5      3];
        [5      4];
        [6      0.5];
        [6      1.5];
        [6      2.5];
        [6      3.5];
        [6      4.5];
        [6      5.5];
        [7      0];
        [7      1];
        [7      2];
        [7      3];
        [7      4];
        [7      5];
        [7      6];
        
        %E
        [10      0];
        [10      1];
        [10      2];
        [10      3];
        [10      4];
        [10      5];
        [10      6];
        [11      0.5];
        [11      1.5];
        [11      2.5];
        [11      3.5];
        [11      4.5];
        [11      5.5];
        [12      0];
        [12      3];
        [12      6];
        [13      0.5];
        [13      2.5];
        [13      5.5];
        [14      0];
        [14      3];
        [14      6];
        
        %L
        [17      0];
        [17      1];
        [17      2];
        [17      3];
        [17      4];
        [17      5];
        [17      6];
        [18      0.5];
        [18      1.5];
        [18      2.5];
        [18      3.5];
        [18      4.5];
        [18      5.5];
        [19      0];
        [19      1];
        [20      0.5];
        [20      1.5];
        
        %L
        [23      0];
        [23      1];
        [23      2];
        [23      3];
        [23      4];
        [23      5];
        [23      6];
        [24      0.5];
        [24      1.5];
        [24      2.5];
        [24      3.5];
        [24      4.5];
        [24      5.5];
        [25      0];
        [25      1];
        [26      0.5];
        [26      1.5];
        
        %O
        [29      2];
        [29      3];
        [29      4];
        [30      0.5];
        [30      1.5];
        [30      4.5];
        [30      5.5];
        [31      0];
        [31      6];
        [32      0.5];
        [32      5.5];
        [33      0];
        [33      6];
        [34      0.5];
        [34      1.5];
        [34      4.5];
        [34      5.5];
        [35      2];
        [35      3];
        [35      4];
        
        ]*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end

function [ Model ] = model_HELLOWORLD( raster_spacing, CCD_pixelsize )

        %HELLO
        Model = [

        %H
        [1      0];
        [1      1];
        [1      2];
        [1      3];
        [1      4];
        [1      5];
        [1      6];
        [2      0.5];
        [2      1.5];
        [2      2.5];
        [2      3.5];
        [2      4.5];
        [2      5.5];
        [3      2];
        [3      3];
        [4      2.5];
        [4      3.5];
        [5      3];
        [5      4];
        [6      0.5];
        [6      1.5];
        [6      2.5];
        [6      3.5];
        [6      4.5];
        [6      5.5];
        [7      0];
        [7      1];
        [7      2];
        [7      3];
        [7      4];
        [7      5];
        [7      6];
        
        %E
        [10      0];
        [10      1];
        [10      2];
        [10      3];
        [10      4];
        [10      5];
        [10      6];
        [11      0.5];
        [11      1.5];
        [11      2.5];
        [11      3.5];
        [11      4.5];
        [11      5.5];
        [12      0];
        [12      3];
        [12      6];
        [13      0.5];
        [13      2.5];
        [13      5.5];
        [14      0];
        [14      3];
        [14      6];
        
        %L
        [17      0];
        [17      1];
        [17      2];
        [17      3];
        [17      4];
        [17      5];
        [17      6];
        [18      0.5];
        [18      1.5];
        [18      2.5];
        [18      3.5];
        [18      4.5];
        [18      5.5];
        [19      0];
        [19      1];
        [20      0.5];
        [20      1.5];
        
        %L
        [23      0];
        [23      1];
        [23      2];
        [23      3];
        [23      4];
        [23      5];
        [23      6];
        [24      0.5];
        [24      1.5];
        [24      2.5];
        [24      3.5];
        [24      4.5];
        [24      5.5];
        [25      0];
        [25      1];
        [26      0.5];
        [26      1.5];
        
        %O
        [29      2];
        [29      3];
        [29      4];
        [30      0.5];
        [30      1.5];
        [30      4.5];
        [30      5.5];
        [31      0];
        [31      6];
        [32      0.5];
        [32      5.5];
        [33      0];
        [33      6];
        [34      0.5];
        [34      1.5];
        [34      4.5];
        [34      5.5];
        [35      2];
        [35      3];
        [35      4];
        
        %W
        [0      -14];
        [0      -13];
        [0      -12];
        [0      -11];
        [0      -10];
        [0      -9];
        [1      -14.5];
        [1      -13.5];
        [2      -15];
        [2      -14];
        [3      -14.5];
        [3      -13.5];
        [3      -12.5];
        [4      -15];
        [4      -14];
        [5      -14.5];
        [5      -13.5];
        [6      -14];
        [6      -13];
        [6      -12];
        [6      -11];
        [6      -10];
        [6      -9];
        
        %O
        [9      -13];
        [9      -12];
        [9      -11];
        [10     -14.5];
        [10     -13.5];
        [10     -10.5];
        [10     -9.5];
        [11     -15];
        [11     -9];
        [12     -14.5];
        [12     -9.5];
        [13     -15];
        [13     -9];
        [14     -14.5];
        [14     -13.5];
        [14     -10.5];
        [14     -9.5];
        [15     -13];
        [15     -12];
        [15     -11];
        
        %R
        [18     -15];
        [18     -14];
        [18     -13];
        [18     -12];
        [18     -11];
        [18     -10];
        [18     -9];
        [19     -14.5];
        [19     -13.5];
        [19     -12.5];
        [19     -11.5];
        [19     -10.5];
        [19     -9.5];
        [20     -13];
        [20     -12];
        [20     -9];
        [21     -14.5];
        [21     -13.5];
        [21     -11.5];
        [21     -9.5];
        [22     -15];
        [22     -14];
        [22     -11];
        [22     -10];
        
        %L
        [25     -15];
        [25     -14];
        [25     -13];
        [25     -12];
        [25     -11];
        [25     -10];
        [25     -9];
        [26     -14.5];
        [26     -13.5];
        [26     -12.5];
        [26     -11.5];
        [26     -10.5];
        [26     -9.5];
        [27     -15];
        [27     -14];
        [28     -14.5];
        [28     -13.5];
        
        %D
        [31     -15];
        [31     -14];
        [31     -13];
        [31     -12];
        [31     -11];
        [31     -10];
        [31     -9];
        [32     -14.5];
        [32     -9.5];
        [33     -15];
        [33     -9];
        [34     -14.5];
        [34     -13.5];
        [34     -10.5];
        [34     -9.5];
        [35     -13];
        [35     -12];
        [35     -11];
        
        ]*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end

function [ Model ] = model_R( raster_spacing, CCD_pixelsize )

        %R
        Model = [
        
        %R
        [18     -15];
        [18     -14];
        [18     -13];
        [18     -12];
        [18     -11];
        [18     -10];
        [18     -9];
        [19     -14.5];
        [19     -13.5];
        [19     -12.5];
        [19     -11.5];
        [19     -10.5];
        [19     -9.5];
        [20     -13];
        [20     -12];
        [20     -9];
        [21     -14.5];
        [21     -13.5];
        [21     -11.5];
        [21     -9.5];
        [22     -15];
        [22     -14];
        [22     -11];
        [22     -10];
              
        ]*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end


function [ Model ] = model_MPI( raster_spacing, CCD_pixelsize )

        Model = [

        [0 0.0];
        [6 0.0];
        [8 0.0];
        [14 0.0];
        [0 1.0];
        [6 1.0];
        [8 1.0];
        [14 1.0];
        [11 1.5];
        [0 2.0];
        [6 2.0];
        [8 2.0];
        [14 2.0];
        [11 2.5];
        [0 3.0];
        [4 3.0];
        [6 3.0];
        [8 3.0];
        [10 3.0];
        [12 3.0];
        [14 3.0];
        [3 3.5];
        [5 3.5];
        [0 4.0];
        [2 4.0];
        [6 4.0];
        [8 4.0];
        [10 4.0];
        [12 4.0];
        [14 4.0];
        [0 5.0];
        [2 5.0];
        [6 5.0];
        [8 5.0];
        [10 5.0];
        [12 5.0];
        [14 5.0];
        [9 5.5];
        [13 5.5];
        [0 6.0];
        [2 6.0];
        [6 6.0];
        [8 6.0];
        [14 6.0];
        [3 6.5];
        [5 6.5];
        [9 6.5];
        [13 6.5];
        [0 7.0];
        [4 7.0];
        [6 7.0];
        [8 7.0];
        [14 7.0];
        
        ]*[-1 0;0 -1].*raster_spacing./CCD_pixelsize;

end

function [ Model ] = model_phantum( raster_spacing, CCD_pixelsize )

        %phantum
        Model = [

        %
        [10     10];
        [10     15];
        [10     20];

        ]*0.2*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end

function [ Model ] = model_phantum_wide( raster_spacing, CCD_pixelsize )

        %phantum
        Model = [

        %
        [10     10];
        [10     20];
        [10     30];

        ]*0.2*[1 0;0 -1].*raster_spacing./CCD_pixelsize;
    
end