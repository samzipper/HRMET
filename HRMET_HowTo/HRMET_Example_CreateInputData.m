%% HRMET_Example_CreateInputData.m
% This script is intended to generate input spatially distributed input
% data to be used in our HRMET example.
%
% In our example, we want to determine what the ET rates were in/around
% Philadelphia PA during the signing of the Declaration of Independence in
% 1776. Note that, by necessity, all input data will be made up.

close all; clear all; clc;    % clean up workspace

%% Create a meshgrid of latitude and longitude. This will be the grid that
% our model runs over, so all other input grids must be the same size.
% For our example problem, we have a square grid, centered on Philadelphia
% (39.95 N, 75.15 W), at 0.00001 degree (approx. 1 m) intervals. The total
% grid size is 41x41.
[long, lat] = meshgrid(75.1498:0.00001:75.1502, 39.9498:0.00001:39.9502);

%% Define input data that does not vary spatially here.
% While some of this input data probably should vary over your model
% domain, in this example, we are treating it as constant.
datetime = 648857.5;   % July 4, 1776, at noon
albSoil  = 0.105;      % Soil albedo
albVeg   = 0.200;      % Vegetation albedo
emissSoil= 0.945;      % Soil emissivity
emissVeg = 0.940;      % Vegetation emissivity
pa       = 101.3;      % Atmospheric pressure [kPa]
LAI      = 2.5;        % Leaf Area Index [m2 m-2] - in the 1770s, an LAI of 2.5 was probably considered quite good.
h        = 1.5;        % Canopy height [m]
SWin     = 700;        % Incoming shortwave radiation [W m-2]
u        = 4.75;       % Wind speed [m s-1]
ea       = 2.25;       % Air vapor pressure [kPa] - this will likely be calculed from temperature & relative humidity in your meteorological dataset

%% Define input data that varies spatially here
% In our example, we will say that canopy temperature and air temperature
% vary spatially.

Tair = 25+randn(size(lat));  % Air temperature [degC] - varies randomly

[X, Y] = meshgrid(linspace(10,15,41), linspace(10,15,41));  % temporary variables to create T input
T = X+Y;                     % Canopy surface temperature [degC] - gradient from 20 in upper left to 30 in lower right
clear X Y

% Generate plots to look at input data, if you want
subplot(2,1,1);
imagesc(Tair); 
colorbar; caxis([20 30]);
title('Air Temperature [C]');

subplot(2,1,2);
imagesc(T);
colorbar; caxis([20 30]);
title('Canopy Surface Temperature [C]');

%% Save workspace
save('HRMET_Example_InputData.mat');
