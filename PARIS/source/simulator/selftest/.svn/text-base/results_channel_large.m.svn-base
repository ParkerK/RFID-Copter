% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_channel_large function
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
filename  = 'results_channel_large';
  
% random test setup bounds (uniformly distributed)
testsettings.ntx  = [  1,    5]; % number of transmitters
testsettings.nrx  = [  1,    5]; % number of receivers
testsettings.dist = [0.1,   10]; % m distance
testsettings.pl   = [1.5,    6]; % path-loss factor
testsettings.c    = [2e8,  4e8]; % m/s speed of light
testsettings.f0   = [1e6, 30e6]; % Hz "center" frequencies for attenuation
testsettings.fs   = [1e9, 30e9]; % Hz sampling frequency

% partitioning
partitioning.names = {'log-dist'}; % == settings.type
partitioning.runs  = [1000];


% *******************************************************************************************************
% complete/check partitioning

if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check partitioning.');
end

% indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% create settings

for i = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)
      % number of transmitters, receivers
      ntx = round(testsettings.ntx(1) - 0.5 + rand*(testsettings.ntx(2) - testsettings.ntx(1)  + 1));
      nrx = round(testsettings.nrx(1) - 0.5 + rand*(testsettings.nrx(2) - testsettings.nrx(1)  + 1));
      % distance matrix (geometry neglected)
      settings{j}.dist = testsettings.dist(1) + rand(nrx,ntx)*(testsettings.dist(2) - testsettings.dist(1) );
      % other settings
      settings{j}.type = partitioning.names{i};
      settings{j}.pl   = testsettings.pl(1)   + rand       *(testsettings.pl(2)   - testsettings.pl(1));
      settings{j}.c    = testsettings.c(1)    + rand       *(testsettings.c(2)    - testsettings.c(1) );
      settings{j}.f0   = testsettings.f0(1)   + rand(1,ntx)*(testsettings.f0(2)   - testsettings.f0(1));
      settings{j}.fs   = testsettings.fs(1)   + rand       *(testsettings.fs(2)   - testsettings.fs(1));
      % round distances to sampling resolution
      settings{j}.dist = round(settings{j}.dist/settings{j}.c*settings{j}.fs)*settings{j}.c/settings{j}.fs;
   end
end


% *******************************************************************************************************
% calculate expected results

for i = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)
      % gain factor (straightforward and thus slightly awkward)
      f0 = repmat(settings{j}.f0, size(settings{j}.dist,1), 1); % replicate f0 for receivers
      results{j}.gain = ( settings{j}.c ./ (4*pi*settings{j}.dist.*f0) ) .^ (settings{j}.pl/2);
      results{j}.gain(results{j}.gain > 1) = 1;
      % delay in samples (not rounded on purpose)
      results{j}.delay_s = settings{j}.dist / settings{j}.c * settings{j}.fs;
   end
end


% *******************************************************************************************************
% complete/check partitioning

% remove first entry from partitioning.indices (not needed in test function)
partitioning.indices(1) = [];

% one last check
if ( sum(partitioning.runs) ~= length(results) ) || ( sum(partitioning.runs) ~= length(settings) )
   criterr('Unexpected length of results or settings array. Check partitioning.');
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
characteristic = 'settings and results for selftest: channel_large.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
