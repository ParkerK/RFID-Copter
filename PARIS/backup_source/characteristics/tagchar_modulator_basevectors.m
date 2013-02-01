% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate data for: tag - modulation (base vectors for assembly and detuning states)
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
disp(sprintf('WARNING: This will overwrite existing base vectors!\n'));
reply = input('Are you REALLY sure y/n [n]?  ', 's');
if ~strcmpi(reply, 'y')
   return
end
disp(' '); clear('reply');

% "switch to here"
cd(fileparts(mfilename('fullpath')));


% *******************************************************************************************************
% define assembly / detuning vectors to be able to quickly classify a characteristic

% assembly (perfect tuning down to massive detuning)
v1 = [50:25:800, 850:50:1000, 1250:250:1750, 2000:1000:5000, 7500, 10000] * 1e-15; % F Cat
v2 = [-98, -95, -90, -80, -67, -50, -33, 0, 50, 100, 200, 500, 1000, 2000, 5000, 10000]; % percent Rat shift
assembly.description = 'Capacity .cat [F], shift of Ohmic part .rats [%] from optimum. Assembly state is index in .v_*.';
assembly.cat  = v1;
assembly.rats = v2;
assembly.v_cat  = sort(repmat(v1, 1, length(v2)));
assembly.v_rats =      repmat(v2, 1, length(v1));

% antenna detuning including tolerances (extremely weak resonance up to massive detuning near metal)
v1 = -0.9 : 0.1 : 1; % boost of resonance (-1: no res., 0: perfect ... 0.3: near water ... 1: near metal)
v2 =  0 : 5e6 : 200e6; % Hz shift of resonance
detuning.description = 'Antenna own-resonance boost .res_en -1..1 (0:none), frequency shift .fshift [Hz]. Detuning state is index in .v_*.';
detuning.res_en = v1;
detuning.fshift = v2;
detuning.v_res_en = sort(repmat(v1, 1, length(v2)));
detuning.v_fshift =      repmat(v2, 1, length(v1));


% *******************************************************************************************************
% save results

% prepare information header
matfilename    = 'tagchar_modulator_adsm'; %#ok<NASGU>
characteristic = 'assembly and detuning state mapping (to settings)'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save
save(fullfile('tag', matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
      'assembly', 'detuning');
