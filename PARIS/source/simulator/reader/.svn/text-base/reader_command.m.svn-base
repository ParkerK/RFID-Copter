% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% reader - create commmand data stream (supported: query, ack)
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
% version = reader_command()
%    Just returns the version number (string).
% data = reader_command(command, settings)
%    Returns the bitstream for COMMAND according to SETTINGS.
%    ATTENTION: SETTINGS not checked for conformity/validity/...!
%
%
% ***** Interface definition *****
% data = reader_command(command, settings)
%    command    string defining the command {'query', 'ack'}
%    settings   struct containing parameters for command (are named like in EPCglobal)
%       .dr        (command: 'query') {8, 64/3}
%       .m         (command: 'query') {1, 2, 4, 8}
%       .trext     (command: 'query') {0, 1}
%       .sel       (command: 'query') {'all', '~sl', 'sl'}
%       .session   (command: 'query') {'s0', 's1', 's2', 's3'}
%       .target    (command: 'query') {'a', 'b'}
%       .q         (command: 'query') {0, 1, 2, ..., 14, 15}
%       .rn16      (command: 'ack')   (hex string without 0x)
%
%    data       binary bitstream for command according to settings
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

function data = reader_command(command, settings)
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
% input parameter checks / prepare input parameters

% number of input parameters
if nargin < 2
   criterr('Not enough input arguments.');
end


% *******************************************************************************************************
% commands

switch lower(command)
   
% query
   case 'query'
      % check contents of settings
      %     prepare required data
      expected.name = 'settings';
      expected.reqfields = {'dr', 'm', 'trext', 'sel', 'session', 'target', 'q'};
      %     check
      errortext = contentcheck(settings, expected);
      %     output
      if ~isempty(errortext)
         err('Incomplete settings\n%s', errortext);
      end
      % create command bitstream
      data = zeros(22,1);
      data(1:4) = [1; 0; 0; 0]; % command code
      switch settings.dr % TRcal divide ratio
         case 8
            data(5) = 0;
         case 64/3
            data(5) = 1;
         otherwise
            err('Unsupported DR "%f"', settings.dr);
      end
      switch settings.m % cycles per symbol
         case 1
            data(6:7) = [0; 0];
         case 2
            data(6:7) = [0; 1];
         case 4
            data(6:7) = [1; 0];
         case 8;
            data(6:7) = [1; 1];
         otherwise
            err('Unsupported M "%f"', settings.m);
      end
      data(8) = settings.trext; % pilot tone
      switch lower(settings.sel) % which tags should respond to query
         case 'all'
            data(9:10) = [0; 0];
         case '~sl'
            data(9:10) = [1; 0];
         case 'sl'
            data(9:10) = [1; 1];
         otherwise
            err('Unsupported sel "%s"', settings.sel);
      end
      switch lower(settings.session) % session for inventory round
         case 's0'
            data(11:12) = [0; 0];
         case 's1'
            data(11:12) = [0; 1];
         case 's2'
            data(11:12) = [1; 0];
         case 's3'
            data(11:12) = [1; 1];
         otherwise
            err('Unsupported target "%s"', settings.target);
      end
      switch lower(settings.target) % inventoried flag A/B
         case 'a'
            data(13) = 0;
         case 'b'
            data(13) = 1;
         otherwise
            err('Unsupported target "%s"', settings.target);
      end
      data(14:17) = bitget(settings.q, 4:-1:1); % q value
      data(18:22) = crc(data(1:17), 5); % CRC-5
      
% ACK
   case 'ack'
      % check contents of settings
      %     prepare required data
      expected.name = 'settings';
      expected.reqfields = {'rn16'};
      %     check
      errortext = contentcheck(settings, expected);
      %     output
      if ~isempty(errortext)
         err('Incomplete settings\n%s', errortext);
      end
      % create command bitstream
      data = zeros(18,1);
      data(1:2)  = [0; 1]; % command code
      data(3:18) = hex2bit(settings.rn16); % RN16
      
   otherwise
      err('Unsupported command "%s"', command);
end


