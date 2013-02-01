% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate data for: tag - power supply (Vdda vs. chip input power)
%
%
%
% ***** Copyright / License / Authors *****
% Copyright 2007, 2008, 2009, 2010, 2011 Daniel Arnitz
%   Signal Processing and Speech Communication Laboratory, Graz University of Technology, Austria
%   NXP Semiconductors Austria GmbH Styria, Gratkorn, Austria
% Copyright 2012 Daniel Arnitz
%   Reynolds Lab, Department of Electrical and Computer Engineering, Duke University, USA
%
% This file is part of the PARIS Simulation Framework.
%
% The PARIS Simulation Framework is free software: you can redistribute it and/or modify it under the
% terms of the GNU General Public License as published by the Free Software Foundation, either
% version 3 of the License, or (at your option) any later version.
%
% The PARIS Simulation Framework is distributed in the hope that it will be useful, but WITHOUT ANY
% WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% See the GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License along with the PARIS Simulation 
% Framework. If not, see <http://www.gnu.org/licenses/>.
%
%
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
%
%
% ***** Todo *****
%
% *******************************************************************************************************


% ******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% one last chance to reconsider
disp(sprintf('WARNING: This will overwrite existing files!\n'));
reply = input('Are you REALLY sure y/n [n]?  ', 's');
if ~strcmpi(reply, 'y')
   return
end
disp(' ');

% paths, etc
%     path to globalinit.m
path = fullfile(fileparts(mfilename('fullpath')), '..');
%     Matlab does not resolve symbolic links; if we're running on Unix, let the system resolve them
if isunix
   [dummy, path] = system(sprintf('readlink -f "%s"', path));
end
%     add this path and call globalinit script
restoredefaultpath; addpath(path); clear('dummy', 'path');
globalinit('init');

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% power drain: RC (parallel) load with known input power and voltage => RC is time constant
settings.fs  =               8192e6; % Hz sampling frequency
settings.tau = 2*pi*333e3*360.6e-12; % s time constant

% data points per dimension
settings.n  = 256;

% output subfolder for plots
settings.plotfolder = 'plots';


% ******************************************************************************************************
% load data
meas_tag_vdda = loadmat('meas_tag_vdda');


% ******************************************************************************************************
% RC lowpass filter coefficients

% use c2d to discretize filter and finally tfdata to get filter coefficients
[rc_b, rc_a] = tfdata(c2d(tf(1, [settings.tau, 1]), 1/settings.fs), 'v');


% ******************************************************************************************************
% vdda(pic) characteristic: interpolate to settings.n points

% add one additional point: vdda = 0 for minimum power
meas_tag_vdda.pic      = [meas_tag_vdda.pic(1) - diff(meas_tag_vdda.pic(1:2)); meas_tag_vdda.pic];
meas_tag_vdda.vdda_v2a = [0; meas_tag_vdda.vdda_v2a];
meas_tag_vdda.vdda_v2b = [0; meas_tag_vdda.vdda_v2b];
meas_tag_vdda.vdda_v2c = [0; meas_tag_vdda.vdda_v2c];
meas_tag_vdda.vdda_avg = [0; meas_tag_vdda.vdda_avg];

% interpolate
pic = interp1(1:length(meas_tag_vdda.pic), meas_tag_vdda.pic, linspace(1,length(meas_tag_vdda.pic),settings.n), 'spline')';
vdda_v2a = interp1(1:length(meas_tag_vdda.pic), meas_tag_vdda.vdda_v2a, linspace(1,length(meas_tag_vdda.pic),settings.n), 'spline')';
vdda_v2b = interp1(1:length(meas_tag_vdda.pic), meas_tag_vdda.vdda_v2b, linspace(1,length(meas_tag_vdda.pic),settings.n), 'spline')';
vdda_v2c = interp1(1:length(meas_tag_vdda.pic), meas_tag_vdda.vdda_v2c, linspace(1,length(meas_tag_vdda.pic),settings.n), 'spline')';
vdda_avg = interp1(1:length(meas_tag_vdda.pic), meas_tag_vdda.vdda_avg, linspace(1,length(meas_tag_vdda.pic),settings.n), 'spline')';


% ******************************************************************************************************
% plot

% verification of interpolation
figure; hold on;
plot(meas_tag_vdda.pic, meas_tag_vdda.vdda_v2a, 'b.');
plot(meas_tag_vdda.pic, meas_tag_vdda.vdda_v2b, 'k.');
plot(meas_tag_vdda.pic, meas_tag_vdda.vdda_v2c, 'g.');
plot(meas_tag_vdda.pic, meas_tag_vdda.vdda_avg, 'r.');
plot(pic, vdda_v2a, 'b-');
plot(pic, vdda_v2b, 'k-');
plot(pic, vdda_v2c, 'g-');
plot(pic, vdda_avg, 'r-');
hold off; grid on;
xlim(xyzlimits(pic)); ylim(xyzlimits(vdda_v2a,vdda_v2b,vdda_v2c,vdda_avg));
setlabels('INTERPOLATION OF vdda(pic)', 'Pic [dBm]', 'Vdda [V]');
setlegend({'Sofie\_V2A','Sofie\_V2B', 'Sofie\_V2C', 'overall average'}, 'SouthEast');
hold off; grid on;
savefigure(gcf, fullfile(settings.plotfolder, 'tagchar_power'), 'png');


% ******************************************************************************************************
% save interpolated result

% prepare everything
%     header
matfilename    = 'tagchar_power'; %#ok<NASGU>
characteristic = 'vdda(pic) for 3 chips plus average, RC filter coeff'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = sprintf('%s (%s)', meas_tag_vdda.matfilename, meas_tag_vdda.createdby); %#ok<NASGU>

% save
save(fullfile('tag',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'pic','vdda_v2a','vdda_v2b','vdda_v2c','vdda_avg','settings');
