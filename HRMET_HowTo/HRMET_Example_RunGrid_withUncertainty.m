%% HRMET_Example_RunGrid.m
% This script is intended to run HRMET many times over a grid while
% randomly sampling inputs from a user-defined Gaussian distribution in
% order to estimate uncertainty.
%
% In our example, we want to determine what the ET rates were in/around
% Philadelphia PA during the signing of the Declaration of Independence in
% 1776. Note that, by necessity, all input data will be made up.

close all; clear all; clc;    % clean up workspace

rng(1)  % set random number seed

%% Load input data
load('HRMET_Example_InputData.mat');   % Load input data (created by HRMET_Example_CreateInputData.m)

% Note that we have 4 inputs that vary spatially: latitude, longitude, air
% temperature, and canopy surface temperature. All other inputs are
% constant over space. 

%% Define uncertainties
% There are several ways to estimate uncertainty for an input, which will
% be shown below. In any case, what we are doing is defining the mean and
% the standard deviation, which will then be used to sample an input when
% we run HRMET.

% If we think a reasonable estimate of error is 10% of the input data, we can
% use that.
u_mean = u;    % Not necessary, but renaming our input u as the mean
u_std = u*0.1; % Uncertainty is 10% of u

% For a gridded input (e.g. canopy surface temperature), a moving window
% approach can be used. This can be shifted to within your 'for' loop that \
% is used to run HRMET for computational efficiency, if you want.
window = 3;            % the size of the moving window you want to use
[Tair_mean, Tair_std] = HRMET_Example_MovingWindow(Tair,window);

% An estimated error can also simply be defined by the user.
T_mean = T;
T_std = ones(size(T))*5;  % uniform std of 5 degrees for all points

%% Loop over grid points and run HRMET - with no error (same as HRMET_Example_RunGrid.m)

ET_noError = NaN(size(lat));          % make an empty grid to hold output

for i = 1:size(lat,1);        % latitude changes vertically
    for j = 1:size(long,2);   % longitude changes horizontally
        
        % Run HRMET here with no error
        ET_noError(i,j) = HRMET_shared(datetime, long(i,j), lat(i,j), Tair(i,j), ...
            SWin, u, ea, pa, LAI, h, T(i,j), albSoil, albVeg, emissSoil, emissVeg);
        % Note that we only have to use (i,j) for inputs that vary spatially
    end
end

%% Incorporate error and run HRMET over grid points
iter = 20;   % number of iterations for each point

ET_mean = NaN(size(lat));    % make empty grid to hold output     
ET_std  = NaN(size(lat));         

for i = 1:size(lat,1);        % latitude changes vertically
    for j = 1:size(long,2);   % longitude changes horizontally
        
        % set up empty vectors to hold output
        ETCalc = NaN(iter,1);
        uCalc = NaN(iter,1);
        TairCalc = NaN(iter,1);
        TCalc = NaN(iter,1);
        
        for k = 1:iter;       % loop through as many iterations as you want
            % randomly generate inputs
            uCalc(k,1) = u_mean + u_std*randn;
            TairCalc(k,1) = Tair_mean(i,j) + Tair_std(i,j)*randn;
            TCalc(k,1) = T_mean(i,j) + T_std(i,j)*randn;
            
            % Run HRMET!
            ETCalc(k,1) = HRMET_shared(datetime, long(i,j), lat(i,j), TairCalc(k,1), ...
                SWin, uCalc(k,1), ea, pa, LAI, h, TCalc(k,1), albSoil, albVeg, emissSoil, emissVeg);
        end
        
        ET_mean(i,j) = mean(ETCalc);  % calculate mean & std of all iterations
        ET_std(i,j) = std(ETCalc);
    end
end

%% Plot output
figure(2)
subplot(2,3,1);
imagesc(Tair); 
colorbar; caxis([20 30]);
title('Input Air Temperature [C]');

subplot(2,3,2);
imagesc(Tair_mean); 
colorbar; caxis([20 30]);
title('Mean Air Temperature [C]');

subplot(2,3,3);
imagesc(Tair_std); 
colorbar;
title('Std Air Temperature [C]');

subplot(2,3,4);
imagesc(ET_noError);
colorbar; caxis([0.3 0.9]);
title('No Error ET Rate [mm hr-1]')

subplot(2,3,5);
imagesc(ET_mean);
colorbar; caxis([0.3 0.9]);
title('Mean ET Rate [mm hr-1]')

subplot(2,3,6);
imagesc(ET_std);
colorbar;
title('Std ET Rate [mm hr-1]')