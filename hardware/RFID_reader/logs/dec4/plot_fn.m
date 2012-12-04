clear all
close all

% radius of the earth in meters
R = 6371000;

filename = '5m.txt';

data = importdata(filename);

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

scatter(x,y)
title('GPS', 'FontSize', 40)
xlabel('Latitude  ','FontSize',20)
ylabel('Longitude  ','FontSize',20)




