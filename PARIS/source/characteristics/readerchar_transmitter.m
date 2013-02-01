% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate data for: reader - transmitter (nonlinear power amplifier)
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


% ******************************************************************************************************
% settings

settings.n          =     256; % number of data points
settings.plotfolder = 'plots'; % subfolder for plots


% ******************************************************************************************************
% amplitude characteristic of nonlinear power amplifier (in this case a bad one)

nl_x = linspace(-1, 1, settings.n)';
nl_y = 2/pi*sin(pi/2*nl_x);


% ******************************************************************************************************
% save interpolated result

% create a plot
figure; 
cplot(nl_x*100, nl_y*100, {'-',1,'blue'});
grid on; xlim(xyzlimits(nl_x*100)); ylim(xyzlimits(nl_y*100)); 
setlabels('AMPLITUDE CHARACTERISTIC OF POWER AMPLIFIER', 'input amplitude [%]', 'output amplitude [%]');
savefigure(gcf, fullfile(settings.plotfolder, 'readerchar_transmitter') , 'png');

% prepare everything
%     header
matfilename    = 'readerchar_transmitter'; %#ok<NASGU>
characteristic = 'symmetrical amplitude characteristic of nonlinear poweramp nl_y(nl_x)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save
save(fullfile('reader',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'nl_x','nl_y', 'settings');

