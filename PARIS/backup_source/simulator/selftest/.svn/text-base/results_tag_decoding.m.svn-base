% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_tag_decoding function
% !!! EXISTING RESULTS WILL BE OVERWRITTEN - USE WITH CARE !!!
%
% Warning: relies on valid data in results_reader_command.mat
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
filename  = 'results_tag_decoding';

% timings (random uniform) [min ; max]
setup.t0        = 1e-6*[    10;     50]; % s
setup.tari      = 1e-6*[  6.25;     25]; % s
setup.x         =      [   0.5;      1]; % times tari
setup.trcal     =      [   1.1;      3]; % times rtcal
setup.pw        =      [ 0.265;  0.525]; % times tari
setup.pwmin     =                  2e-6; % s
setup.delimiter = 1e-6*[11.875; 13.125]; % s

% tag clock (sampling rate)
setup.fclk = 1.92e6; % Hz


% *******************************************************************************************************
% get data sequences from results_reader_command
rcommand_data = loadmat('results_reader_command', '     ');
settings = rcommand_data.settings;
results  = rcommand_data.results;

% add command information and data streams to settings
for i = 1 : length(settings.query)
   settings.query{i}.command = 'query';
   settings.query{i}.data    = results.query(i, :);
end
for i = 1 : length(settings.ack)
   settings.ack{i}.command = 'ack';
   settings.ack{i}.data    = results.ack(i, :);
end

% reformat
results = [settings.query; settings.ack];
clear settings;


% *******************************************************************************************************
% create timings in [samples] = results for decoding

for i = 1 : length(results)
   % generate random timings
   results{i}.t0        = setup.t0(1)        + rand * ( setup.t0(2)        - setup.t0(1)        );
   results{i}.tari      = setup.tari(1)      + rand * ( setup.tari(2)      - setup.tari(1)      );
   results{i}.x         = setup.x(1)         + rand * ( setup.x(2)         - setup.x(1)         );
   results{i}.trcal     = setup.trcal(1)     + rand * ( setup.trcal(2)     - setup.trcal(1)     );
   results{i}.pw        = setup.pw(1)        + rand * ( setup.pw(2)        - setup.pw(1)        );
   results{i}.delimiter = setup.delimiter(1) + rand * ( setup.delimiter(2) - setup.delimiter(1) );
   
   % times tari
   results{i}.x  = results{i}.x  * results{i}.tari;
   results{i}.pw = results{i}.pw * results{i}.tari;
   
   % restrictions
   if results{i}.pw < setup.pwmin
      results{i}.pw = setup.pwmin;
   end
      
   % us -> samples
   results{i}.t0_s        = round( results{i}.t0        * setup.fclk );
   results{i}.tari_s      = round( results{i}.tari      * setup.fclk );
   results{i}.x_s         = round( results{i}.x         * setup.fclk );
   results{i}.pw_s        = round( results{i}.pw        * setup.fclk );
   results{i}.delimiter_s = round( results{i}.delimiter * setup.fclk );
   
   % rtcal_s and trcal_s are calculated in samples to make sure that rounding (...nonlinear) is no problem
   results{i}.rtcal_s = 2 * results{i}.tari_s + results{i}.x_s;
   results{i}.trcal_s = round( results{i}.rtcal_s * results{i}.trcal );
   
   % samples -> s
   results{i}.t0        = results{i}.t0_s        / setup.fclk;
   results{i}.tari      = results{i}.tari_s      / setup.fclk;
   results{i}.x         = results{i}.x_s         / setup.fclk;
   results{i}.rtcal     = results{i}.rtcal_s     / setup.fclk;
   results{i}.trcal     = results{i}.trcal_s     / setup.fclk;
   results{i}.pw        = results{i}.pw_s        / setup.fclk;
   results{i}.delimiter = results{i}.delimiter_s / setup.fclk;
   
end


% *******************************************************************************************************
% create "demodulated" signals to decode

signals = cell(size(results));

for i = 1 : length(results)
   % patterns @ fclk
   patterns.t0 = ones(results{i}.t0_s, 1);
   patterns.delimiter = zeros(results{i}.delimiter_s, 1);
   patterns.rtcal = [ones(results{i}.rtcal_s - results{i}.pw_s, 1); zeros(results{i}.pw_s ,1)];
   patterns.trcal = [ones(results{i}.trcal_s - results{i}.pw_s, 1); zeros(results{i}.pw_s ,1)];
   patterns.data0 = [ones(results{i}.tari_s  - results{i}.pw_s, 1); zeros(results{i}.pw_s ,1)];
   patterns.data1 = [ones(results{i}.tari_s + results{i}.x_s - results{i}.pw_s, 1); zeros(results{i}.pw_s ,1)];
   
   % create leadin
   if strcmpi(results{i}.command, 'query')
      leadin = [patterns.delimiter; patterns.data0; patterns.rtcal; patterns.trcal]; % preamble
   else
      leadin = [patterns.delimiter; patterns.data0; patterns.rtcal]; % framesync
      % TRcal not needed in results any more (just confusing)
      results{i}.trcal   = NaN;
      results{i}.trcal_s = NaN;
   end
   
   % create data stream
   data = [];
   for j = 1 : length(results{i}.data)
      if results{i}.data(j) == 0
         data = [data; patterns.data0];
      else
         data = [data; patterns.data1];
      end
   end
   
   % assemble
   signals{i} = [patterns.t0; leadin; data; patterns.t0];
   
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
characteristic = 'settings and results for selftest: tag_decoding.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = sprintf('%s (%s)', rcommand_data.matfilename, rcommand_data.createdby); %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'signals', 'results');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
