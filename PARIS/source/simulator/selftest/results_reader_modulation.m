% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_reader_modulation function
% !!! EXISTING RESULTS WILL BE OVERWRITTEN - USE WITH CARE !!!
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
% = add SSB-ASK after overshoot problem is fixed
%
% *******************************************************************************************************


% *******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% initialize global stuff, print only warnings and errors
globalinit('exceptions');

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% warn if "untrustworthy" functions are used
reply = input('Relies on reader_modulation to sanitize settings. Results have to be hand-checked! OK y/n [n]?  ', 's');
if ~strcmpi(reply, 'y')
   return
end


% *******************************************************************************************************
% test setup

% output directory and filename
directory = 'results';
filename  = 'results_reader_modulation';

% (necessary) carrier setup
stdsettings.oscillator.fcenter      =  910e6; % Hz
stdsettings.oscillator.oversampling =     16; % times carrier frequency
% => sampling (HF) frequency
fs = 2*stdsettings.oscillator.fcenter*stdsettings.oscillator.oversampling; % Hz

% standard (initial) setup
stdsettings.modulation.forcesettings =        false; % if true, no consistency/conformity checks are done
stdsettings.modulation.modulation    =    'DSB-ASK'; % {'DSB-ASK', 'SSB-ASK', 'PR-ASK'}
stdsettings.modulation.leadin        =   'preamble'; % {'preamble', 'frame-sync' or 'framesync'}
stdsettings.modulation.tari          =      6.25e-6; % s duration of data-0 (preferred: {6.25, 12.5, 25 us})
stdsettings.modulation.lf            =        640e3; % Hz backscatter link frequency {40...640 kHz}
stdsettings.modulation.dr            =         64/3; % divide ratio DR {64/3, 8}
stdsettings.modulation.delimiter     =        13e-6; % s (initial) delimiter duration; -1: random
stdsettings.modulation.x             =          0.5; % times tari {0.5...1 ; -1 means random} (data-1=tari+x)
stdsettings.modulation.pw            =          0.4; % times tari {0.265 ... 0.525 ; -1 means random}
stdsettings.modulation.moddepth      =           90; % percent {80...100 ; -1 means random}
stdsettings.modulation.trf           =         0.32; % times tari {0...0.33 ; -1 means random}
stdsettings.modulation.t0            =         2e-6; % s time delay (unmodulated carrier before and after modulation)
stdsettings.modulation.fs            =           fs; % Hz

% data (random sequence length taken from this vector)
stdsettings.data.length = [2, 4, 6, 8]; % has to be even
      
% nonrandom test setups (will also lead to consistency problems)
testsettings.modulation = {'DSB-ASK', 'PR-ASK'};
testsettings.leadin     = {'preamble', 'frame-sync'};
testsettings.tari       = [6, 8.111, 12.5, 30]*1e-6; % s
testsettings.lf         = [20, 333, 700]*1e3; % Hz
testsettings.dr         = [64/3, 8];
testsettings.delimiter  = linspace(12.5e-6*0.9, 12.5e-6*1.1, 3); % s
testsettings.x          = linspace(0.4, 1.3, 3); % times tari
testsettings.pw         = linspace(0.1, 0.5, 3); % times tari
testsettings.moddepth   = linspace(70, 100, 3); % percent
testsettings.trf        = linspace(0.05, 0.4, 3); % % times tari (problem for timing detection if <<)
testsettings.t0         = linspace(2e-6, 10e-6, 3); % s (problem for timing detection if <<)

% partitioning (has to match test setup)
partitioning.names = {'inconsistent settings', 'forced inconsistent settings',... 
   'sweep modulation', 'sweep leadin', 'sweep tari',...
   'sweep lf', 'sweep dr', 'sweep delimiter', 'sweep x', 'sweep pw', 'sweep moddepth',...
   'sweep trf', 'sweep t0', 'random'};
partitioning.runs  = [1, 1, 2, 2, 4, 3, 2, 3, 3, 3, 3, 3, 3, 32];

% carrier setup
oscsettings.fcenter      =   910e6; % Hz
oscsettings.fstddev      =    10e3; % Hz
oscsettings.astddev      =    1e-7; % amplitude instability (variance)
oscsettings.snr          =      80; % dB
oscsettings.fs           = 16*2*oscsettings.fcenter; % Hz (16x oversampling)


% *******************************************************************************************************
% create settings: sweeps

% standard setup, violations
settings{1} = stdsettings.modulation;
settings{1}.tari          =     6e-6; % s
settings{1}.lf            = 169.31e3; % Hz
settings{1}.delimiter     =  17.2e-6; % s

% force settings, violations
settings{2} = stdsettings.modulation;
settings{2}.tari          =     6e-6; % s
settings{2}.lf            = 169.31e3; % Hz
settings{2}.delimiter     =  17.2e-6; % s
settings{2}.forcesettings = true;

% index in settings
index = 3;

% sweep modulation
for i = 1 : length(testsettings.modulation);
   settings{index} = stdsettings.modulation;
   settings{index}.modulation = testsettings.modulation{i};
   index = index + 1;   
end

% sweep leadin
for i = 1 : length(testsettings.leadin);
   settings{index} = stdsettings.modulation;
   settings{index}.leadin = testsettings.leadin{i};
   index = index + 1;   
end

% sweep tari
for i = 1 : length(testsettings.tari);
   settings{index} = stdsettings.modulation;
   settings{index}.tari = testsettings.tari(i);
   index = index + 1;   
end

% sweep lf
for i = 1 : length(testsettings.lf);
   settings{index} = stdsettings.modulation;
   settings{index}.lf = testsettings.lf(i);
   index = index + 1;   
end

% sweep dr
for i = 1 : length(testsettings.dr);
   settings{index} = stdsettings.modulation;
   settings{index}.dr = testsettings.dr(i);
   index = index + 1;   
end

% sweep delimiter
for i = 1 : length(testsettings.delimiter);
   settings{index} = stdsettings.modulation;
   settings{index}.delimiter = testsettings.delimiter(i);
   index = index + 1;   
end

% sweep x
for i = 1 : length(testsettings.x);
   settings{index} = stdsettings.modulation;
   settings{index}.x = testsettings.x(i);
   index = index + 1;   
end

% sweep pw
for i = 1 : length(testsettings.pw);
   settings{index} = stdsettings.modulation;
   settings{index}.pw = testsettings.pw(i);
   index = index + 1;   
end

% sweep moddepth
for i = 1 : length(testsettings.moddepth);
   settings{index} = stdsettings.modulation;
   settings{index}.moddepth = testsettings.moddepth(i);
   index = index + 1;   
end

% sweep trf
for i = 1 : length(testsettings.trf);
   settings{index} = stdsettings.modulation;
   settings{index}.trf = testsettings.trf(i);
   index = index + 1;   
end

% sweep t0
for i = 1 : length(testsettings.t0);
   settings{index} = stdsettings.modulation;
   settings{index}.t0 = testsettings.t0(i);
   index = index + 1;   
end


% *******************************************************************************************************
% create settings: random

for i = index : index + partitioning.runs(end) - 1
   settings{i} = stdsettings.modulation;
   settings{i}.delimiter  = -1;
   settings{i}.x          = -1;
   settings{i}.pw         = -1;
   settings{i}.moddepth   = -1;
   settings{i}.trf        = -1;
end


% *******************************************************************************************************
% create settings: data with random length (only even: PR-ASK)

for i = 1 : length(settings)
   % length
   n = stdsettings.data.length( 1 + round(rand*(length(stdsettings.data.length)-1)) );
   % create data stream (make sure there is at least one 0 and one 1)
   data = [];
   while sum(data) == 0 || sum(data) == n
      data = round( rand(n, 1) );
   end
   settings{i}.data = data;
end


% *******************************************************************************************************
% create results ("result" = completed and sanitized setup) and complete/check partitioning

% max( results{i}.length ) for carrier creation (speeds up the process)
oscsettings.length = 0;

% create
for i = 1 : length(settings)
   % initialize RNGs (used by reader_modulation and reader_oscillator)
   rand('twister', i);
   randn('state', i);
   % determine needed carrier length, complete and sanitize settings
   results{i} = reader_modulation(settings{i}.data, settings{i});
   % add data to results
   results{i}.data = settings{i}.data;
   % find longest length_s (for carrier creation)
   if oscsettings.length <= results{i}.length
      oscsettings.length = results{i}.length;
   end
end

% partitioning
if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check partitioning.');
end
if sum(partitioning.runs) ~= length(results)
   criterr('Unexpected length of results array. Check partitioning.');
end
partitioning.indices = cumsum(partitioning.runs); % index of end-of-block


% *******************************************************************************************************
% save results

% is there already a .mat-file with this name?
search = dir(strcat(fullfile(directory, filename),'*'));
if ~isempty( search )
   % ask before overwriting
   reply = input('This will overwrite an existing file. Are you REALLY sure y/n [n]?  ', 's');
   if ~strcmpi(reply, 'y')
      return
   end
end

% prepare information header
matfilename    = filename; %#ok<NASGU>
characteristic = 'settings and results for selftest: reader_modulation.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning', 'oscsettings');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
