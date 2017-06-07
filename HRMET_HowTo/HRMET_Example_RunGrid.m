%% HRMET_Example_RunGrid.m
% This script is intended to run HRMET over a grid and plot the output.
%
% In our example, we want to determine what the ET rates were in/around
% Philadelphia PA during the signing of the Declaration of Independence in
% 1776. Note that, by necessity, all input data will be made up.

close all; clear all; clc;    % clean up workspace

%% Load input data
load('HRMET_Example_InputData.mat');   % Load input data (created by HRMET_Example_CreateInputData.m)

% Note that we have 4 inputs that vary spatially: latitude, longitude, air
% temperature, and canopy surface temperature. All other inputs are
% constant over space. 

%% Loop over grid points and run HRMET

ET = NaN(size(lat));          % make an empty grid to hold output

for i = 1:size(lat,1);        % latitude changes vertically
    for j = 1:size(long,2);   % longitude changes horizontally
        
        % Run HRMET here
        ET(i,j) = HRMET_shared(datetime, long(i,j), lat(i,j), Tair(i,j), ...
            SWin, u, ea, pa, LAI, h, T(i,j), albSoil, albVeg, emissSoil, emissVeg);
        % Note that we only have to use (i,j) for inputs that vary spatially
    end
end

%% Plot output
subplot(3,1,1);
imagesc(Tair); 
colorbar; caxis([20 30]);
title('Air Temperature [C]');

subplot(3,1,2);
imagesc(T);
colorbar; caxis([20 30]);
title('Canopy Surface Temperature [C]');

subplot(3,1,3);
imagesc(ET);
colorbar; caxis([0.3 0.9]);
title('Instantaneous ET Rate [mm hr-1]')

% As you can see, the ET rate is strongly controlled by the temperature
% gradient from the canopy surface to the atmosphere. In the lower-right,
% where the canopy surface temperature is higher than the air temperature,
% lots of energy is partitioned to the sensible heat flux. In the contrast,
% in the upper-left, the air temperature is generally warmer than the
% canopy surface, and more energy is partitioned to latent heat flux.