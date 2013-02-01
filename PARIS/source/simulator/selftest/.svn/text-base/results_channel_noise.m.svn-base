% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_channel_noise function
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
filename  = 'results_channel_noise';
   
% random test setup bounds (uniformly distributed)

testsettings.size = [5e5,  1e6; 1,   5]; % [samples length of test vector; number of receivers]
testsettings.n0   = [  0, -100; 1, 1e5]; % [dBm; Hz] single-sided noise density [dBm, per ? Hz] @ receiver
testsettings.p_id =         0.1; % probability for identical setups for all receivers
testsettings.fs   = [1e9, 30e9]; % Hz sampling frequency
testsettings.frxs = [1e6, 30e6]; % Hz receiver (tag,reader) sampling frequency (for noise density -> power)

% partitioning
partitioning.names = {'off', 'wgn'}; % == settings.type
partitioning.runs  = [   10,   100];


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
      % general
      settings{j}.size = round([testsettings.size(1,1) - 0.5 + rand*(testsettings.size(1,2) - testsettings.size(1,1) + 1),...
                                testsettings.size(2,1) - 0.5 + rand*(testsettings.size(2,2) - testsettings.size(2,1) + 1)]);
      settings{j}.type = partitioning.names{i};
      settings{j}.fs   = testsettings.fs(1)   + rand*(testsettings.fs(2)   - testsettings.fs(1)  );
      % identical variance for all receivers?
      if round(rand-0.5+testsettings.p_id)
         % yes: identical n0 and frxs
         settings{j}.n0   = [testsettings.n0(1,1) + rand*(testsettings.n0(1,2) - testsettings.n0(1,1)),...
                             testsettings.n0(2,1) + rand*(testsettings.n0(2,2) - testsettings.n0(2,1))];
         settings{j}.frxs = testsettings.frxs(1) + rand*(testsettings.frxs(2) - testsettings.frxs(1));
      else
         % no: n0 and frxs i.i.d. (per receiver)
         settings{j}.n0   = [testsettings.n0(1,1) + rand(settings{j}.size(2),1)*(testsettings.n0(1,2) - testsettings.n0(1,1)),...
                             testsettings.n0(2,1) + rand(settings{j}.size(2),1)*(testsettings.n0(2,2) - testsettings.n0(2,1))];
         settings{j}.frxs = testsettings.frxs(1) + rand(settings{j}.size(2),1)*(testsettings.frxs(2) - testsettings.frxs(1));
      end
   end
end


% *******************************************************************************************************
% calculate expected results

for i = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)
      % noise variance
      if strcmpi(settings{j}.type, 'wgn') % settings{j}.n0 is [dBm; Hz]
         results{j}.var = 10.^(settings{j}.n0(:,1)/10)*1e-3./settings{j}.n0(:,2) * settings{j}.fs./settings{j}.frxs;
      else
         results{j}.var = zeros(1, settings{j}.size(2));
      end
      %     make sure size is correct in case of identical setups
      if settings{j}.size(2) > 1 && length(results{j}.var) == 1
         results{j}.var = repmat(results{j}.var, 1, settings{j}.size(2));
      end
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
characteristic = 'settings and results for selftest: channel_noise.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
