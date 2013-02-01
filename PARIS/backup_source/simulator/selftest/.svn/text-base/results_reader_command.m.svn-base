% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_reader_command function
% !!! EXISTING RESULTS WILL BE OVERWRITTEN - USE WITH CARE !!!
%
% WARNING: Results have to be checked manually (relies on reader_command)!
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

% warn if "untrustworthy" functions are used
reply = input('Relies on reader_command. Results have to be hand-checked! OK y/n [n]?  ', 's');
if ~strcmpi(reply, 'y')
   return
end


% *******************************************************************************************************
% test setup

% output directory and filename
directory = 'results';
filename  = 'results_reader_command';

% query command
query.runs     = 1000;
query.dr       = [8, 64/3];
query.m        = [1, 2, 4, 8];
query.trext    = [0, 1];
query.sel      = {'all', '~sl', 'sl'};
query.session  = {'s0', 's1', 's2', 's3'};
query.target   = {'a', 'b'};
query.q        = 0:1:15;

% ack command
ack.runs = 1000;
ack.rn16 = round(linspace(0,2^16-1,97)); % (97 is prime)


% *******************************************************************************************************
% initialize arrays

% query command
settings.query = cell(query.runs, 1);
results.query  = zeros(query.runs, 22);

% ack command
settings.ack  = cell(ack.runs, 1);
results.ack   = zeros(ack.runs, 18);


% *******************************************************************************************************
% create settings: 27 sweeps for query command

index = 1;

% initial setup for sweeps (the following loops depend on initial setting to index 1)
sweepset.dr       = query.dr(1);
sweepset.m        = query.m(1);
sweepset.trext    = query.trext(1);
sweepset.sel      = query.sel{1};
sweepset.session  = query.session{1};
sweepset.target   = query.target{1};
sweepset.q        = query.q(1);

% sweep dr
for i = 1 : length(query.dr);
   sweepset.dr = query.dr(i);
   settings.query{index} = sweepset;
   index = index + 1;
end

% sweep m
for i = 2 : length(query.m);
   sweepset.m = query.m(i);
   settings.query{index} = sweepset;
   index = index + 1;
end

% sweep trext
for i = 2 : length(query.trext);
   sweepset.trext = query.trext(i);
   settings.query{index} = sweepset;
   index = index + 1;
end

% sweep sel
for i = 2 : length(query.sel);
   sweepset.sel = query.sel{i};
   settings.query{index} = sweepset;
   index = index + 1;
end

% sweep session
for i = 2 : length(query.session);
   sweepset.session = query.session{i};
   settings.query{index} = sweepset;
   index = index + 1;
end

% sweep target
for i = 2 : length(query.target);
   sweepset.target = query.target{i};
   settings.query{index} = sweepset;
   index = index + 1;
end

% sweep q
for i = 2 : length(query.q);
   sweepset.q = query.q(i);
   settings.query{index} = sweepset;
   index = index + 1;
end


% *******************************************************************************************************
% create settings: query.runs-27 random query commands

for i = index : query.runs
   % generate random sequence for setup selection
   randseq.dr      = 1 + round(rand * (length(query.dr)-1) );
   randseq.m       = 1 + round(rand * (length(query.m)-1) );
   randseq.trext   = 1 + round(rand * (length(query.trext)-1) );
   randseq.sel     = 1 + round(rand * (length(query.sel)-1) );
   randseq.session = 1 + round(rand * (length(query.session)-1) );
   randseq.target  = 1 + round(rand * (length(query.target)-1) );
   randseq.q       = 1 + round(rand * (length(query.q)-1) );
   % select settings 
   settings.query{i}.dr      = query.dr(randseq.dr);
   settings.query{i}.m       = query.m(randseq.m);
   settings.query{i}.trext   = query.trext(randseq.trext);
   settings.query{i}.sel     = query.sel{randseq.sel};
   settings.query{i}.session = query.session{randseq.session};
   settings.query{i}.target  = query.target{randseq.target};
   settings.query{i}.q       = query.q(randseq.q); % in case .q is not 0:1:15 (does not really matter)
end


% *******************************************************************************************************
% create settings: ack command

% 97 runs linear
for i = 1 : length(ack.rn16)
   settings.ack{i}.rn16 = dec2hex(ack.rn16(i), 4);
end

% all the rest: random
for i = length(ack.rn16)+1 : ack.runs
   settings.ack{i}.rn16 = dec2hex(round(rand*(2^16-1)), 4);
end


% *******************************************************************************************************
% create results (hand-check!)

% query command
for i = 1 : query.runs
   results.query(i,:) = reader_command('query', settings.query{i});
end

% ack command
for i = 1 : ack.runs
   results.ack(i,:) = reader_command('ack', settings.ack{i});
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
characteristic = 'settings and results for selftest: reader_command.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
              'settings', 'results');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
