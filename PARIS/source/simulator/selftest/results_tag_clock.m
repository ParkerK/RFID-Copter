% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_tag_clock function
% !!! EXISTING RESULTS WILL BE OVERWRITTEN - USE WITH CARE !!!
%
% WARNING: Results have to be checked manually (results created by tag_clock)!
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
% tests setup

% output directory and filename
directory = 'results';
filename  = 'results_tag_clock';

% runs per mode
tests.runs        = 100;

% borders for setups
tests.fs          = 1e6*[500, 10000]; % [min, max] sampling frequency in Hz
tests.fcenter     = 1e6*[0.5,    10]; % [min, max] center frequency in Hz
tests.fsigma      =     [  0,    10]; % [min, max] sampling time variation (sigma) in percent
tests.phi0        =     [  0,   360]; % [min, max] degree phase delay
tests.length_fs   =     [1e7,   1e8]; % [min, max] max(clock) for mode = 'fs';
tests.length_fclk =     [1e3,   1e5]; % [min, max] length(clock) for mode = 'fclk';

% minimum values (everything smaller will be set to zero to avoid problems with rel. error)
tests.fsigma_min  = 0.2; % percent
tests.phi0_min    =   5; % degree


% *******************************************************************************************************
% create settings: mode = 'fs'

% counters for special events
phi0_rand   = 0;
phi0_zero   = 0;
fsigma_zero = 0;

% create
for i = 1 : tests.runs
   % create random setup
   settings.mode_fs(i).mode      = 'fs';
   settings.mode_fs(i).fs        = tests.fs(1)        + rand * (tests.fs(2)        - tests.fs(1)     );
   settings.mode_fs(i).fcenter   = tests.fcenter(1)   + rand * (tests.fcenter(2)   - tests.fcenter(1));
   settings.mode_fs(i).fsigma    = tests.fsigma(1)    + rand * (tests.fsigma(2)    - tests.fsigma(1) );
   settings.mode_fs(i).phi0      = tests.phi0(1)      + rand * (tests.phi0(2)      - tests.phi0(1)   );
   settings.mode_fs(i).length    = round(tests.length_fs(1) + rand * (tests.length_fs(2) - tests.length_fs(1)));
     
   % chance for random setup of t0: 25%
   if rand > 0.75 % random
      rand('twister', i);  % initialize RNG for reproducible results
      settings.mode_fs(i).phi0 = -1; % set to random
      settings.mode_fs(i).phi0_res = mod(rand*360, 360); % expected result
      phi0_rand = phi0_rand + 1;
   else
      settings.mode_fs(i).phi0_res = settings.mode_fs(i).phi0; 
   end
   
   % minimum values (avoid problems with rel. error)
   if settings.mode_fs(i).phi0_res < tests.phi0_min
      settings.mode_fs(i).phi0     = 0;
      settings.mode_fs(i).phi0_res = 0;
      phi0_zero = phi0_zero + 1;
   end
   if settings.mode_fs(i).fsigma < tests.fsigma_min
      settings.mode_fs(i).fsigma = 0;
      fsigma_zero = fsigma_zero + 1;
   end
end

% output
disp(sprintf('phi0: %i random, %i zero  ;  fsigma: %i zero   (everything should be nonzero here!)', phi0_rand, phi0_zero, fsigma_zero));


% *******************************************************************************************************
% create settings: mode = 'fs'

% counters for special events
phi0_rand   = 0;
phi0_zero   = 0;
fsigma_zero = 0;

% create
for i = 1 : tests.runs
   % create random setup
   settings.mode_fclk(i).mode      = 'fclk';
   settings.mode_fclk(i).fs        = tests.fs(1)          + rand * (tests.fs(2)          - tests.fs(1)      );
   settings.mode_fclk(i).fcenter   = tests.fcenter(1)     + rand * (tests.fcenter(2)     - tests.fcenter(1) );
   settings.mode_fclk(i).fsigma    = tests.fsigma(1)      + rand * (tests.fsigma(2)      - tests.fsigma(1)  );
   settings.mode_fclk(i).phi0      = tests.phi0(1)        + rand * (tests.phi0(2)        - tests.phi0(1)    );
   settings.mode_fclk(i).length    = round(tests.length_fclk(1) + rand * (tests.length_fclk(2) - tests.length_fclk(1)));
   
   % chance for random setup of t0: 25%
   if rand > 0.75 % random
      rand('twister', 12345+i);  % initialize RNG for reproducible results
      settings.mode_fclk(i).phi0 = -1; % set to random
      settings.mode_fclk(i).phi0_res = mod(rand*360, 360); % expected result
      phi0_rand = phi0_rand + 1;
   else
      settings.mode_fclk(i).phi0_res = settings.mode_fclk(i).phi0; 
   end
   
   % minimum values (avoidproblems with rel. error)
   if settings.mode_fclk(i).phi0_res < tests.phi0_min
      settings.mode_fclk(i).phi0     = 0;
      settings.mode_fclk(i).phi0_res = 0;
      phi0_zero = phi0_zero + 1;
   end
   if settings.mode_fclk(i).fsigma < tests.fsigma_min
      settings.mode_fclk(i).fsigma = 0;
      fsigma_zero = fsigma_zero + 1;
   end
end

% output
disp(sprintf('phi0: %i random, %i zero  ;  fsigma: %i zero   (everything should be nonzero here!)', phi0_rand, phi0_zero, fsigma_zero));


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
characteristic = 'settings and results for selftest: tag_clock.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
   
