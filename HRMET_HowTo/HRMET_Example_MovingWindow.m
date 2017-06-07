function [out_mean, out_std] = HRMET_Example_MovingWindow(in, window)
%% HRMET_Example_MovingWindow.m
% This script is intended to pass a moving window of size 'window' over a
% dataset 'in' and return the mean and standard deviation for each out in
% 'out'.

shift = (window-1)/2;  % this is distance the edge of your window is from the center

out_mean = NaN(size(in)); % empty grid to hold output
out_std =  NaN(size(in)); 
for i = 1+shift:size(in,1)-shift;   % shift the start of your moving window away from the edges, because it cannot be calculated there
    for j = 1+shift:size(in,2)-shift;
        out_mean(i,j)= mean2(in(i-shift:i+shift, j-shift:j+shift));
        out_std(i,j) = std2(in(i-shift:i+shift, j-shift:j+shift));
    end
end
