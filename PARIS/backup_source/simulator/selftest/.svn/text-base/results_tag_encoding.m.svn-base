% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_tag_encoding function
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
% initialization and settings
clear; close all; clc;
version = 'beta 3.0';

% initialize global stuff, print only warnings and errors (assumes that globalinit.m is part of path)
globalinit('exceptions');

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% output directory and filename
directory = 'results';
filename  = 'results_tag_encoding';


% *******************************************************************************************************
% expected results

% FM0 encoding (rising edge at beginning => init at s3 or s4)
fm0.seq00 = -[1,0,1,0]*2+1;
fm0.seq01 = -[1,0,1,1]*2+1;
fm0.seq10 = -[1,1,0,1]*2+1;
fm0.seq11 = -[1,1,0,0]*2+1;

% FM0 preamble
fm0.preamble_trext0 = [1,1,0,1,0,0,1,0,0,0,1,1]*2-1;
fm0.preamble_trext1 = [repmat([1,0]*2-1,1,12), fm0.preamble_trext0];

% FM0 End-of-Signaling
fm0.eos0 = [1,1]*2-1; % last "halfbit" was low
fm0.eos1 = [0,0]*2-1; % last "halfbit" was high

% Miller encoding, M=2 (rising edge at beginning => init at s1 or s2)
miller2.seq000 = [1,0,1,0, 0,1,0,1, 1,0,1,0]*2-1;
miller2.seq001 = [1,0,1,0, 0,1,0,1, 0,1,1,0]*2-1;
miller2.seq010 = [1,0,1,0, 1,0,0,1, 0,1,0,1]*2-1;
miller2.seq011 = [1,0,1,0, 1,0,0,1, 0,1,1,0]*2-1;
miller2.seq100 = [1,0,0,1, 0,1,0,1, 1,0,1,0]*2-1;
miller2.seq101 = [1,0,0,1, 0,1,0,1, 0,1,1,0]*2-1;
miller2.seq110 = [1,0,0,1, 0,1,1,0, 1,0,1,0]*2-1;
miller2.seq111 = [1,0,0,1, 0,1,1,0, 1,0,0,1]*2-1;

% Miller preamble, M=2
miller2.preamble_trext0 = [repmat([1,0]*2-1,1, 4*2), miller2.seq010, -miller2.seq111];
miller2.preamble_trext1 = [repmat([1,0]*2-1,1,16*2), miller2.seq010, -miller2.seq111];

% Miller End-of-Signaling, M=2
miller2.eos = [1,0,0,1]*2-1;

% Miller encoding, M=4 (rising edge at beginning => init at s1 or s2)
miller4.seq000 = [1,0,1,0,1,0,1,0, 0,1,0,1,0,1,0,1, 1,0,1,0,1,0,1,0]*2-1;
miller4.seq001 = [1,0,1,0,1,0,1,0, 0,1,0,1,0,1,0,1, 0,1,0,1,1,0,1,0]*2-1;
miller4.seq010 = [1,0,1,0,1,0,1,0, 1,0,1,0,0,1,0,1, 0,1,0,1,0,1,0,1]*2-1;
miller4.seq011 = [1,0,1,0,1,0,1,0, 1,0,1,0,0,1,0,1, 0,1,0,1,1,0,1,0]*2-1;
miller4.seq100 = [1,0,1,0,0,1,0,1, 0,1,0,1,0,1,0,1, 1,0,1,0,1,0,1,0]*2-1;
miller4.seq101 = [1,0,1,0,0,1,0,1, 0,1,0,1,0,1,0,1, 0,1,0,1,1,0,1,0]*2-1;
miller4.seq110 = [1,0,1,0,0,1,0,1, 0,1,0,1,1,0,1,0, 1,0,1,0,1,0,1,0]*2-1;
miller4.seq111 = [1,0,1,0,0,1,0,1, 0,1,0,1,1,0,1,0, 1,0,1,0,0,1,0,1]*2-1;

% Miller preamble, M=4
miller4.preamble_trext0 = [repmat([1,0]*2-1,1, 4*4), miller4.seq010, -miller4.seq111];
miller4.preamble_trext1 = [repmat([1,0]*2-1,1,16*4), miller4.seq010, -miller4.seq111];

% Miller End-of-Signaling, M=4
miller4.eos = [1,0,1,0, 0,1,0,1]*2-1;

% Miller encoding, M=8 (rising edge at beginning => init at s1 or s2)
miller8.seq000 = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0, 0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]*2-1;
miller8.seq001 = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0, 0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0]*2-1;
miller8.seq010 = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0, 1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1]*2-1;
miller8.seq011 = [1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0, 1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0]*2-1;
miller8.seq100 = [1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]*2-1;
miller8.seq101 = [1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0]*2-1;
miller8.seq110 = [1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0, 1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0]*2-1;
miller8.seq111 = [1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1, 0,1,0,1,0,1,0,1,1,0,1,0,1,0,1,0, 1,0,1,0,1,0,1,0,0,1,0,1,0,1,0,1]*2-1;

% Miller preamble, M=8
miller8.preamble_trext0 = [repmat([1,0]*2-1,1, 4*8), miller8.seq010, -miller8.seq111];
miller8.preamble_trext1 = [repmat([1,0]*2-1,1,16*8), miller8.seq010, -miller8.seq111];

% Miller End-of-Signaling, M=8
miller8.eos = [1,0,1,0, 1,0,1,0, 0,1,0,1, 0,1,0,1]*2-1;


% *******************************************************************************************************
% assemble test patterns - FM0 encoding

% no pilot (TRext=0)
results.fm0.trext0 =...
[fm0.preamble_trext0, fm0.seq00, fm0.eos1;...
 fm0.preamble_trext0, fm0.seq01, fm0.eos0;...
 fm0.preamble_trext0, fm0.seq10, fm0.eos0;...
 fm0.preamble_trext0, fm0.seq11, fm0.eos1]';

% pilot (TRext=1)
results.fm0.trext1 =...
[fm0.preamble_trext1, fm0.seq00, fm0.eos1;...
 fm0.preamble_trext1, fm0.seq01, fm0.eos0;...
 fm0.preamble_trext1, fm0.seq10, fm0.eos0;...
 fm0.preamble_trext1, fm0.seq11, fm0.eos1]';


% *******************************************************************************************************
% assemble test patterns - Miller encoding

% M=2, no pilot (TRext=0)
results.miller2.trext0 =...
[miller2.preamble_trext0, miller2.seq000, miller2.eos;...
 miller2.preamble_trext0, miller2.seq001, miller2.eos;...
 miller2.preamble_trext0, miller2.seq010, -miller2.eos;...
 miller2.preamble_trext0, miller2.seq011, miller2.eos;...
 miller2.preamble_trext0, miller2.seq100, miller2.eos;...
 miller2.preamble_trext0, miller2.seq101, miller2.eos;...
 miller2.preamble_trext0, miller2.seq110, miller2.eos;...
 miller2.preamble_trext0, miller2.seq111, -miller2.eos]';

% M=2, pilot (TRext=1)
results.miller2.trext1 =...
[miller2.preamble_trext1, miller2.seq000, miller2.eos;...
 miller2.preamble_trext1, miller2.seq001, miller2.eos;...
 miller2.preamble_trext1, miller2.seq010, -miller2.eos;...
 miller2.preamble_trext1, miller2.seq011, miller2.eos;...
 miller2.preamble_trext1, miller2.seq100, miller2.eos;...
 miller2.preamble_trext1, miller2.seq101, miller2.eos;...
 miller2.preamble_trext1, miller2.seq110, miller2.eos;...
 miller2.preamble_trext1, miller2.seq111, -miller2.eos]';

% M=4, no pilot (TRext=0)
results.miller4.trext0 =...
[miller4.preamble_trext0, miller4.seq000, miller4.eos;...
 miller4.preamble_trext0, miller4.seq001, miller4.eos;...
 miller4.preamble_trext0, miller4.seq010, -miller4.eos;...
 miller4.preamble_trext0, miller4.seq011, miller4.eos;...
 miller4.preamble_trext0, miller4.seq100, miller4.eos;...
 miller4.preamble_trext0, miller4.seq101, miller4.eos;...
 miller4.preamble_trext0, miller4.seq110, miller4.eos;...
 miller4.preamble_trext0, miller4.seq111, -miller4.eos]';

% M=4, pilot (TRext=1)
results.miller4.trext1 =...
[miller4.preamble_trext1, miller4.seq000, miller4.eos;...
 miller4.preamble_trext1, miller4.seq001, miller4.eos;...
 miller4.preamble_trext1, miller4.seq010, -miller4.eos;...
 miller4.preamble_trext1, miller4.seq011, miller4.eos;...
 miller4.preamble_trext1, miller4.seq100, miller4.eos;...
 miller4.preamble_trext1, miller4.seq101, miller4.eos;...
 miller4.preamble_trext1, miller4.seq110, miller4.eos;...
 miller4.preamble_trext1, miller4.seq111, -miller4.eos]';

% M=8, no pilot (TRext=0)
results.miller8.trext0 =...
[miller8.preamble_trext0, miller8.seq000, miller8.eos;...
 miller8.preamble_trext0, miller8.seq001, miller8.eos;...
 miller8.preamble_trext0, miller8.seq010, -miller8.eos;...
 miller8.preamble_trext0, miller8.seq011, miller8.eos;...
 miller8.preamble_trext0, miller8.seq100, miller8.eos;...
 miller8.preamble_trext0, miller8.seq101, miller8.eos;...
 miller8.preamble_trext0, miller8.seq110, miller8.eos;...
 miller8.preamble_trext0, miller8.seq111, -miller8.eos]';

% M=8, pilot (TRext=1)
results.miller8.trext1 =...
[miller8.preamble_trext1, miller8.seq000, miller8.eos;...
 miller8.preamble_trext1, miller8.seq001, miller8.eos;...
 miller8.preamble_trext1, miller8.seq010, -miller8.eos;...
 miller8.preamble_trext1, miller8.seq011, miller8.eos;...
 miller8.preamble_trext1, miller8.seq100, miller8.eos;...
 miller8.preamble_trext1, miller8.seq101, miller8.eos;...
 miller8.preamble_trext1, miller8.seq110, miller8.eos;...
 miller8.preamble_trext1, miller8.seq111, -miller8.eos]';


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
characteristic = 'settings and results for selftest: tag_encoding.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'results');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
