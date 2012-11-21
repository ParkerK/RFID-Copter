clear all
close all

% radius of the earth in meters
R = 6371000;

filename = 'modified/full_spiral_end_on_tag.txt';
%filename = 'modified/first_run.txt';
%filename = 'modified/fast_fly_by.txt';
%filename = 'modified/medium_fly_by.txt';
%filename = 'modified/slow_fly_by.txt';
%filename = 'modified/first_run.txt';

data = importdata(filename);
data(:,[1,3]) = data(:,[3,1]);

% In Durham:
% Lat_Degreess approx 35
% Lat_Minutes approx 59
% Lat_Seconds approx 54.927
% Long_Degrees approx -78
% Long_Minutes approx 56
% Long_Seconds approx 14.8344

lon = data(:,1) ./ 10000000.0;
lat = data(:,2) ./ 10000000.0;

x = R .* cosd(lat) .* cosd(lon);
y = R .* cosd(lat) .* sind(lon);

%filtered_data = imfilter(avg_rssi_matrix, fspecial('average',[2 2]));

avg_rssi_plot = stem3(x,y,data(:,3))
title('RSSI per GPS Coordinate  ', 'FontSize', 40)
zlabel('RSSI','FontSize',20)
xlabel('Latitude  ','FontSize',20)
ylabel('Longitude  ','FontSize',20)




