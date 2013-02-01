% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_channel_small function
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

% maximum matrix size (prevent heavy swapping)
maxram = 6; % gigabyte (approximate) maximum RAM usage during test

% output directory and filename
directory = 'results';
filename  = 'results_channel_small';

% nonrandom setup
stdsettings.seed = round(mod(sum(1e7*clock), 2^32-1) * rand);

% random setup bounds (uniformly distributed)
% ... 
testsettings.maxrays  = [  250,  750]; % maximum number of paths
testsettings.maxiter  = [  250, 1000]; % maximum number of iterations for rms delay spread                  
testsettings.k        = [  -10,   40]; % dB Ricean K-factor
testsettings.trms     = [ 2e-9, 2e-8]; % s RMS delay spread
testsettings.bw       = [  1e7, 1e8 ]; % Hz bandwidth for channel (with oversampling)
testsettings.eps_k    = [ 1e-5, 1e-3]; % relative tolerable error in K-factor (linear)  ...
testsettings.eps_trms = [ 1e-3, 3e-3]; % relative tolerable error in RMS delay spread  ...
testsettings.fres     = [  5e5,  1e7]; % samples/Hz initial frequency resolution
testsettings.fs       = [  1e9, 30e9]; % Hz sampling frequency
testsettings.p_nraysf =          0.03; % probability for forced number of paths (det. LOS only)

% other test settings
testsettings.ensembles = 1e5; % number of independent ensembles for 'rand' averaging

% partitioning
% ... place seed tests RIGHT AFTER non-seed-tests
partitioning.names = {'off', 'deterministic NLOS', 'deterministic NLOS + seed',...
   'random NLOS', 'random NLOS + seed'};
partitioning.runs  = [10, 1000, 100, 10, 5];


% *******************************************************************************************************
% complete/check partitioning

if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check partitioning.');
end

% indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% create settings

% maximum ram usage
ramusage = 0;

for i = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)
      % on/off, number of paths
      settings{j}.on = ~strcmpi(partitioning.names{i}, 'off');
      % length/iteration settings
      settings{j}.maxrays = round(testsettings.maxrays(1) + rand*(testsettings.maxrays(2) - testsettings.maxrays(1) ));
      settings{j}.maxiter = round(testsettings.maxiter(1) + rand*(testsettings.maxiter(2) - testsettings.maxiter(1) ));
      % forced number of paths (can overwrite maxrays)?
      if strcmpi(partitioning.names{i}, 'deterministic NLOS') && rand < testsettings.p_nraysf
         settings{j}.nrays_f = round( rand * 2 * settings{j}.maxrays );
      end
      % other settings
      settings{j}.k        = testsettings.k(1)        + rand*(testsettings.k(2)        - testsettings.k(1)       );
      settings{j}.trms     = testsettings.trms(1)     + rand*(testsettings.trms(2)     - testsettings.trms(1)      );
      settings{j}.bw       = testsettings.bw(1)       + rand*(testsettings.bw(2)       - testsettings.bw(1)        );
      settings{j}.fres     = testsettings.fres(1)     + rand*(testsettings.fres(2)     - testsettings.fres(1)      );
      settings{j}.fs       = testsettings.fs(1)       + rand*(testsettings.fs(2)       - testsettings.fs(1)        );
      settings{j}.eps_k    = testsettings.eps_k(1)    + rand*(testsettings.eps_k(2)    - testsettings.eps_k(1)   );
      settings{j}.eps_trms = testsettings.eps_trms(1) + rand*(testsettings.eps_trms(2) - testsettings.eps_trms(1));
      % number of independent ensembles
      %    ... settings.ensembles is optional for channel_small, but not for test_channel_small
      if any(strcmpi(partitioning.names{i}, {'random NLOS', 'random NLOS + seed'}))
         settings{j}.det       =                  false; % random ...
         settings{j}.ensembles = testsettings.ensembles; %  ... averaging necessary
      else
         settings{j}.det       = true; % deterministic ...
         settings{j}.ensembles =    1; %  ... 1 run is enough
      end
      % deterministic seed tests: identical settings
      if any(strcmpi(partitioning.names{i}, {'deterministic NLOS + seed', 'random NLOS + seed'}))
         settings{j}.seed      = stdsettings.seed;
         settings{j}.ensembles = settings{j-1}.ensembles; % (just to make sure)
         settings{j}.maxrays   = settings{j-1}.maxrays;
         settings{j}.maxiter   = settings{j-1}.maxiter;
         settings{j}.k         = settings{j-1}.k;
         settings{j}.trms      = settings{j-1}.trms;
         settings{j}.bw        = settings{j-1}.bw;
         settings{j}.fres      = settings{j-1}.fres;
         settings{j}.fs        = settings{j-1}.fs;
         settings{j}.eps_k     = settings{j-1}.eps_k;
         settings{j}.eps_trms  = settings{j-1}.eps_trms;
      else
         settings{j}.seed = NaN; % random
      end
      
      % record maximum RAM usage (approximately)
      if ramusage < 4 * settings{j}.maxrays * settings{j}.ensembles  * 8 / 1024^3 % GByte
         ramusage = 4 * settings{j}.maxrays * settings{j}.ensembles  * 8 / 1024^3;
      end
   end
end

% prevent heavy swapping (only for last test)
if  ramusage > maxram
   err('This would need more than %.1f GB of RAM. Reduce maxrays or ensembles.', ramusage);
else
   disp(sprintf('\nPlease note: RAM usage will be about %.1f GB.\n', ramusage));
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
characteristic = 'settings and results for selftest: channel_large.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));

