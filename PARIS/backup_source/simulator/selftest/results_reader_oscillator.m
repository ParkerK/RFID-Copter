% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_reader_oscillator function
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
%
% *******************************************************************************************************


% *******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% initialize global stuff, print only warnings and errors (assumes that globalinit.m is part of path)
globalinit('exceptions');

% "switch to here"
cd(fileparts(mfilename('fullpath')));


% *******************************************************************************************************
% test setup

% output directory and filename
directory = 'results';
filename  = 'results_reader_oscillator';
 
% random test setup bounds (uniformly distributed)
testsettings.fcenter = [  1e8,  1e9]; % Hz
testsettings.fstddev = [  1e0,  1e6]; % Hz (min>>eps)
testsettings.astddev = [ 1e-6, 1e-2]; % V  (min>>eps)
testsettings.snr     = [  100,   10]; % dB
testsettings.length  = [ 1e-5, 5e-5]; % s
testsettings.osr     = [  1.5,    5]; % oversampling rate for fs
testsettings.mode    = {'sin', 'cos', 'exp'}; % carrier type

% partitioning
partitioning.names = {'clean carrier', 'amplitude instability', 'frequency instability', 'additive noise'};
partitioning.runs  = [100, 100, 100, 100];


% *******************************************************************************************************
% complete/check partitioning

if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check partitioning.');
end

% indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% create settings

for j = 1 : length(partitioning.names)
   for i = 1 + partitioning.indices(j) : partitioning.indices(j+1)
      % define clean carrier
      settings{i}.mode    = testsettings.mode{randi([1, length(testsettings.mode)])};
      settings{i}.fcenter = testsettings.fcenter(1) + rand*(diff(testsettings.fcenter));
      settings{i}.osr     = testsettings.osr(1) + rand*(diff(testsettings.osr));
      settings{i}.fs      = 2 * settings{i}.osr * settings{i}.fcenter;
      settings{i}.length  = testsettings.length(1)  + rand*(diff(testsettings.length ));
      settings{i}.fstddev =   0;
      settings{i}.astddev =   0;
      settings{i}.snr     = Inf;
      % add instabilities as required by the current test
      switch lower(partitioning.names{j})
         case 'clean carrier' % just to allow for "otherwise"
         case 'amplitude instability'
            settings{i}.astddev = testsettings.astddev(1) + rand*(diff(testsettings.astddev));
         case 'frequency instability'
            settings{i}.fstddev = testsettings.fstddev(1) + rand*(diff(testsettings.fstddev));
            % modify fs to account for higher frequencies due to unfiltered random frequency component
            settings{i}.fs      = 2 * settings{i}.osr * (settings{i}.fcenter + 5*settings{i}.fstddev);
         case 'additive noise'
            settings{i}.snr     = testsettings.snr(1)     + rand*(diff(testsettings.snr    ));
         otherwise
            error('Unsupported test "%s".', lower(partitioning.names{j}));
      end
   end
end


% *******************************************************************************************************
% create (missing) results

for i = 1 : length(settings)
   % amplitude
   results{i}.ampl = sqrt(2);
   % phase (degree)
   switch lower(settings{i}.mode)
      case 'sin'
         results{i}.phase_re =  0;
         results{i}.phase_im =  0;
      case 'cos'
         results{i}.phase_re = 90;
         results{i}.phase_im = 90;
      case 'exp'
         results{i}.phase_re = 90;
         results{i}.phase_im =  0;
      otherwise
         error('Unsupported settings{%i}.mode="%s".', i, lower(settings{i}.mode));
   end
end


% *******************************************************************************************************
% complete/check partitioning

% remove first entry from partitioning.indices (not needed in test function)
partitioning.indices(1) = [];

% one last check
if ( sum(partitioning.runs) ~= length(settings) )
   criterr('Unexpected length of settings array. Check partitioning.');
end


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
characteristic = 'settings and results for selftest: reader_oscillator.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles', ...
               'settings','results', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
