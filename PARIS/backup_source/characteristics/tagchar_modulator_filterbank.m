% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate data for: tag - modulation (filterbank prototype filter coefficients)
%
% WARNING: Settings are not checked for sanity!
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

% *******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% paths
path = fullfile(fileparts(mfilename('fullpath')), '..');
globalinit('silent');

% one last chance to reconsider
disp(sprintf('WARNING: This will overwrite existing coefficients!\n'));
reply = input('Are you REALLY sure y/n [n]?  ', 's');
if ~strcmpi(reply, 'y')
   return
end
disp(' '); clear('reply');

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% settings
settings.nfilt     =    32; % filter coefficients per channel filter (FIR)
settings.nfft      =  4096; % channels (set to settings.nf for filterbank with frequency warping)


% *******************************************************************************************************
% prepare filter coefficients for filterbank channel filters ("prototype")
disp(sprintf('Preparing filter coefficients for filterbank channel filters ("prototypes")...'));

coeff = npr_coeff(settings.nfft, settings.nfilt);


% *******************************************************************************************************
% save results
disp(sprintf('\nSaving results...'));

% prepare information header
matfilename    = 'tagchar_modulator_fb'; %#ok<NASGU>
characteristic = 'filterbank channel filter coefficients'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save
save(fullfile('tag', matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
      'coeff','settings');
