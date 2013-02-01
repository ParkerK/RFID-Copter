% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: tag - decoding
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
% ***** Behavior *****
% version = test_tag_decoding()
%    Just returns the version number (string).
% sumoferrors = test_tag_decoding(sumoferrors)
%    Tests the function tag_decoding, displays the results in the command window and returns the 
%    overall amount of errors.
%
%
% ***** Interface definition *****
% function sumoferrors = test_reader_modulation(sumoferrors);
%    sumoferrors    sum of errors found not including this fcn call
%
%    sumoferrors    sum of errors found including this fcn call (+1 for each erroneous tested
%                   functionality)
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


function sumoferrors = test_tag_decoding(sumoferrors)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   sumoferrors = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% settings for tag_decoding
settings.fclk = 1.92e6; % tag clock in Hz

% partitioning
partitioning.names = {'query command', 'ack command'};
partitioning.runs  = [1000, 1000];


% *******************************************************************************************************
% run tests and evaluate

% output
disp('   = tag_decoding **');
disp('     checking command, data including CRC, timings (+/-1%)');

% load expected results
data = loadmat('results_tag_decoding', '     ');
signals  = data.signals;
results  = data.results;

% partitioning
if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check internal setup of this fcn.');
end
if sum(partitioning.runs) ~= length(results)
   criterr('Unexpected length of results array. Please check this fcn and results_reader_modulation.mat');
end
partitioning.indices = cumsum(partitioning.runs); % index of end-of-block

% run tests
for i = 1 : length(signals)
   
   % partitioning and output
   %     find index for partition struct
   for j = 1 : length(partitioning.indices)
      if i <= partitioning.indices(j)
         index = j;
         break;
      end
   end
   
   % display once for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
   end
   
   % run tag_decoding
   [data, linkinfo] = tag_decoding(signals{i}, settings);
   
   % reset errors once for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      errors.command = 0;
      errors.timing  = 0;
      errors.data    = 0;
      errors.crc5_ok = 0;
      errortext      = ''; % makes it easier to identify the source of a timing error
   end
        
   % check timings
   %     tari
   [errors.timing, errortext] = check_result(i, linkinfo.tari_s, results{i}.tari_s, errors.timing,...
      errortext, 'tari_s', 'relerr', 0.01);
   [errors.timing, errortext] = check_result(i, linkinfo.tari,   results{i}.tari,   errors.timing,...
      errortext, 'tari',   'relerr', 0.01);
   %     rtcal
   [errors.timing, errortext] = check_result(i, linkinfo.rtcal_s, results{i}.rtcal_s, errors.timing,...
      errortext, 'rtcal_s', 'relerr', 0.01);
   [errors.timing, errortext] = check_result(i, linkinfo.rtcal,   results{i}.rtcal,   errors.timing,...
      errortext, 'rtcal',   'relerr', 0.01);
   %     trcal (only for query command ... preamble)
   if strcmpi(linkinfo.cmd, 'query')
      [errors.timing, errortext] = check_result(i, linkinfo.trcal_s, results{i}.trcal_s, errors.timing,...
         errortext, 'trcal_s', 'relerr', 0.01);
      [errors.timing, errortext] = check_result(i, linkinfo.trcal,   results{i}.trcal  , errors.timing,...
         errortext, 'trcal',   'relerr', 0.01);
   end
   %     threshold between data-0 and data-1
   [errors.timing, errortext] = check_result(i, linkinfo.threshold_s, results{i}.tari_s+results{i}.x_s/2, errors.timing,...
      errortext, 'threshold', 'relerr', 0.01);
   
   % check decoded command, data and crc5
   if strcmpi(linkinfo.cmd, results{i}.command) % type of command ok?
      % command information
      switch lower(linkinfo.cmd)
         case 'query'
            % command info
            if    linkinfo.dr    ~= results{i}.dr    || ...
                  linkinfo.m     ~= results{i}.m     || ...
                  linkinfo.trext ~= results{i}.trext || ...
                  linkinfo.q     ~= results{i}.q     || ...
                  ~strcmpi(linkinfo.sel, results{i}.sel) || ...
                  ~strcmpi(linkinfo.session, results{i}.session) || ...
                  ~strcmpi(linkinfo.target,  results{i}.target)
               errors.command = errors.command + 1;
            end
            % crc5_ok (has to be ok)
            errors.crc5_ok = errors.crc5_ok + (~linkinfo.crc5_ok);
         case 'ack'
            % just RN16
            errors.command = errors.command + (~strcmpi(linkinfo.rn16, results{i}.rn16));
            
         otherwise
            criterr('Unsupported command in results and linkinfo!');
      end
      % decoded data (including crc-5 if present)
      errors.data = check_result(i, data, results{i}.data, errors.data, '', '', 'equal');
   else
      errors.command = errors.command + 1;
   end
   
   % output (if current loop is last test in block)
   if i == partitioning.indices(index)
      if (errors.command + errors.timing + errors.data + errors.crc5_ok) == 0
         disp('         ... passed');
      else
         disp(sprintf('         ... ERRORS: %i command, %i timing, %i data, %i crc5_ok',...
            errors.command, errors.timing, errors.data, errors.crc5_ok));
         sumoferrors = sumoferrors + 1;  
      end
   end

end

