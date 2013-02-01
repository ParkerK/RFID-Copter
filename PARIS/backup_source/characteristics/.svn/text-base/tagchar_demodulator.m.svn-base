% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate data for: tag - demodulation (AGC characteristics)
%
% Only characteristic for venv is used because venv-vpeak = 50mV except for the saturation in vpeak.
% This saturation can easily be implemented by a simple if structure (no lookup table needed).
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
addpath(path); clear('dummy', 'path');
globalinit('init');

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% RC lowpass
settings.rc_fs  = [2e9, 1e8, 32e9]; % [min, delta, max] Hz sampling frequency
settings.rc_tau = 169.12e-9; % s time constant (RC)

% AGC characteristic: data points per dimension
settings.n_vrf     = 96;
settings.n_vdda    = 64;
settings.vdda_snip = 42; % cut interpolation including/not including vdda=1.25V

% other
settings.plotfolder    = 'plots'; % subfolder for plots


% ******************************************************************************************************
% RC lowpass filter coefficients

% frequency vector
rc_fs = [settings.rc_fs(1) : settings.rc_fs(2) : settings.rc_fs(3)]';
settings.n_rc = length(rc_fs);

% use c2d to discretize filter and finally tfdata to get filter coefficients
for i = 1 : settings.n_rc
   [rc_b(i,:), rc_a(i,:)] = tfdata(c2d(tf(1, [settings.rc_tau, 1]), 1/rc_fs(i)), 'v'); %#ok<*SAGROW>
end


% ******************************************************************************************************
% AGC characteristic: interpolate to settings.n_vrf x settings.n_vdda
% (twice for good linear behavior with vdda for vdda > 1.5 and smooth for vdda < 1.5)

% load data
meas_tag_agc = loadmat('meas_tag_agc');

% interpolate without vdda=1.25V (approx linear with vdda)
%     axes
vrf_i  = linspace(meas_tag_agc.vrf(1),  meas_tag_agc.vrf(end),  settings.n_vrf);
vdda_i = linspace(meas_tag_agc.vdda(2), meas_tag_agc.vdda(end), settings.vdda_snip);
%     grid for interp2 (Vrf is line)
[X, Y]     = meshgrid(meas_tag_agc.vdda(2:end), meas_tag_agc.vrf);
[X_i, Y_i] = meshgrid(vdda_i, vrf_i);
%     interpolate
venv_i1 = interp2(X,Y, meas_tag_agc.venv(:, 2:end), X_i,Y_i, 'spline');

% interpolate including vdda=1.25V (creates a "wave" in vdda direction)
%     axes
vrf_i  = linspace(meas_tag_agc.vrf(1),  meas_tag_agc.vrf(end),  settings.n_vrf);
vdda_i = linspace(meas_tag_agc.vdda(1), meas_tag_agc.vdda(end), settings.n_vdda);
%     grid for interp2
[X, Y]     = meshgrid(meas_tag_agc.vdda, meas_tag_agc.vrf);
[X_i, Y_i] = meshgrid(vdda_i, vrf_i);
%     interpolate
venv_i = interp2(X,Y, meas_tag_agc.venv, X_i,Y_i, 'spline');

% combine both (remove wave for vdda>=1.5V)
venv_i(:, settings.n_vdda-settings.vdda_snip+1:end) = venv_i1;


% ******************************************************************************************************
% plots

% AGC: verification of interpolation
figure; hold on;
surface(meas_tag_agc.vrf, meas_tag_agc.vdda, meas_tag_agc.venv', 'FaceColor', 'none');
surface(vrf_i, vdda_i, venv_i', 'EdgeColor', 'none');
hold off; grid on; view([-60, 15]);
setlabels('INTERPOLATION OF venv(vrf, vdda)', 'vrf [V] (peak)', 'vdda [V]', 'venv [V] (peak)');
setlegend({'original data', 'interpolated'}, 'NorthWest');
savefigure(gcf, fullfile(settings.plotfolder, 'tagchar_demodulator_venv'), 'png');

% AGC: venv-vpeak
figure; hold on;
surface(meas_tag_agc.vrf, meas_tag_agc.vdda, meas_tag_agc.venv', 'FaceColor', 'none');
surface(meas_tag_agc.vrf, meas_tag_agc.vdda, meas_tag_agc.vpeak', 'EdgeColor', 'none');
hold off; grid on; view([-60, 15]);
setlabels('VENV - VPEAK (not interpolated)', 'vrf [V] (peak)', 'vdda [V]', 'venv [V] (peak)');
setlegend({'venv', 'vpeak'}, 'NorthWest');
savefigure(gcf, fullfile(settings.plotfolder, 'tagchar_demodulator_venv-vpeak'), 'png');


% ******************************************************************************************************
% save interpolated result

% prepare everything
%     header
matfilename    = 'tagchar_demodulator'; %#ok<NASGU>
characteristic = 'venv(vrf, vdda), RC filter coeff rc_a(rc_fs), rc_b(rc_fs)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = sprintf('%s (%s)', meas_tag_agc.matfilename, meas_tag_agc.createdby); %#ok<NASGU>
%     data
vrf = vrf_i(:);
vdda = vdda_i(:);
venv = venv_i;

% save
save(fullfile('tag',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'vrf','vdda','venv', 'rc_fs','rc_a','rc_b', 'settings');

