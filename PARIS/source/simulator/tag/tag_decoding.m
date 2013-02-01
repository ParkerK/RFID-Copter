% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% tag - decoding
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
% version = tag_clock()
%    Just returns the version number (string).
% [data, linkinfo] = tag_decoding(signal, settings)
%    Decodes the demodulated signal (SIGNAL) according to settings and returns the derived information
%    (command, crc-check, timings) in LINKINFO along with the decoded bitstream in DATA. LINKINFO.TRCAL
%    will be NaN if the decoded command is not a query command. Returns empty DATA and a dummy LINKINFO 
%   (NaNs for all timings and cmd='') if nothing was detected in decoded signal.
%    
%
%
% ***** Interface definition *****
% function [data, linkinfo] = tag_decoding(signal, settings)
%    signal
%    settings
%       .fclk      Tag clock (center) frequency in Hz
%
%    data
%    linkinfo   decoded command (including all info fields) and timing information
%       .tari_s (.tari)     Tari in samples (s)
%       .rtcal_s (.rtcal)   RTcal in samples (s)
%       .trcal_s (.trcal)   TRcal in samples (s)
%       .threshold_s        time between two rising edges > threshold: data-1, else data-0
%       .cmd                decoded command
%       -> All other fields depend on the decoded command. See help reader_command for details.
%    
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


function [data, linkinfo] = tag_decoding(signal, settings)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   data = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks

% number of input parameters
if nargin < 2
   criterr('Not enough input arguments.');
end

% length of vectors
if isempty(signal)
   err('Length of input signal is zero.')
end

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'fclk'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end


% *******************************************************************************************************
% scan input data vector and get timings

% get to first falling edge (start of delimiter)
startindex = 1;
while ~(signal(startindex) == 1 && signal(startindex+1) == 0)
   startindex = startindex + 1;
end

% detect rising edges
rising_edges = [];
for i = startindex : length(signal) - 1
   if signal(i) == 0 && signal(i+1) == 1
      rising_edges = [rising_edges; i];
   end
end

% delays between rising edges (in samples)
delays = diff(rising_edges);

% safety check
if length(delays) < 6 % tari plus rtcal plus 4 bit for largest command code (query) to avoid errors
   warn('Did not detect any supported command in demodulated signal.');
   data = [];
   linkinfo = struct('tari_s',NaN, 'tari', NaN, 'rtcal_s', NaN, 'rtcal',NaN, 'trcal_s', NaN,...
      'trcal', NaN, 'threshold_s', NaN, 'cmd','');
   return
end

% obtain necessary data (in periods clk frequency)
linkinfo.tari_s  = delays(1);
linkinfo.rtcal_s = delays(2);
if delays(3) > delays(2)
   linkinfo.trcal_s = delays(3);
   datastart = 4;
else
   linkinfo.trcal_s = NaN;
   datastart = 3;
end

% the same stuff in s
linkinfo.tari  = linkinfo.tari_s  / settings.fclk;
linkinfo.rtcal = linkinfo.rtcal_s / settings.fclk;
linkinfo.trcal = linkinfo.trcal_s / settings.fclk;


% *******************************************************************************************************
% simple threshold detection

% threshold between length of data-0 (tari) and data-1 (rtcal-tari)
linkinfo.threshold_s = linkinfo.rtcal_s / 2;
data = double(delays(datastart:end) > linkinfo.threshold_s);
data = data(:);


% *******************************************************************************************************
% decode command

% query command
if all(data(1:4) == [1;0;0;0]) && (length(data) >= 22)
   linkinfo.cmd     = 'query'; % command
   linkinfo.dr      = bit2opt(data(5), {8, 64/3}); % TRcal divide ratio DR
   linkinfo.m       = bit2opt(data(6:7), {1, 2, 4, 8}); % cycles per symbol M
   linkinfo.trext   = data(8); % pilot tone
   linkinfo.sel     = bit2opt(data(9:10), {'All', 'All', '~SL', 'SL'}); % which tags should respond to query
   linkinfo.session = bit2opt(data(11:12), {'S0', 'S1', 'S2', 'S3'}); % session for inventory round
   linkinfo.target  = bit2opt(data(13), {'A', 'B'}); % inventoried flag A/B
   linkinfo.q       = bit2opt(data(14:17)); % q value
   linkinfo.crc5    = data(18:22); % CRC-5
   linkinfo.crc5_ok = ( sum(crc(data, 5)) == 0 ); % CRC-5 check

% ack command
elseif all(data(1:2) == [0; 1]) && (length(data) >= 18)
   linkinfo.cmd  = 'ack';
   linkinfo.rn16 = dec2hex(bit2opt(data(3:18)), 4);

% unsupported command => ignore and issue a warning
else
   linkinfo.cmd = '';
   warn('Truncated signal, unsupported command, or TRcal <= RTcal (not allowed).');
end

