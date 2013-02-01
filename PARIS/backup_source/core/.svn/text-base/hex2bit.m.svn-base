% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% hex string to bit array (ignoring spaces)
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
% bitarray = hex2bit(hexstring)
%    Returns the binary representation of the hexadecimal number (without any prefix or suffix like '0x')
%    in a bit array. Ignores spaces (allows grouping).
%    e.g. hex2bit('0A 8f') returns [0;0;0;0;1;0;1;0;1;0;0;0;1;1;1;1]
%         hex2bit('')      returns an empty array
%         hex2bit('  ')    returns an empty array
%
%
% ***** Interface definition *****
% function bitarray = hex2bit(hexstring)
%    hexstring    hexadecimal number (no prefix/suffix like '0x')
%
%    bitarray     column vector containing the binary representation of hexstring
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

function bitarray = hex2bit(hexstring)

% *******************************************************************************************************
% input parameter checks

% number of input parameters
if nargin < 1
   error('Not enough input arguments.');
end

% empty hexstring
if isempty(hexstring)
   bitarray = [];
   return
end

% check for invalid characters
validcharsfound = textscan(hexstring, '%[ 0123456789abcdefABCDEF]', 'Whitespace', '');
validcharsfound = validcharsfound{:};
if isempty(validcharsfound) || ( cellfun(@length, validcharsfound) ~= length(hexstring) )
   error('hexstring contains non-hex characters');
end


% *******************************************************************************************************
% hex -> bit array

% remove whitespaces
blocks = textscan(hexstring, '%[0123456789abcdefABCDEF]', 'Whitespace', '', 'Delimiter', ' ', 'MultipleDelimsAsOne', 1);
blocks = blocks{:};
hexstring = '';
for i = 1 : length(blocks)
    hexstring = [hexstring, blocks{i}];
end

% check validity of hex string
if sum(isstrprop(hexstring, 'xdigit'))  ~= length(hexstring)
   error('hexstring contains non-hex characters.')
end

% create bitarray
bitarray = zeros(4*length(hexstring), 1);

% translate hex per hex, ignore spaces
for i = 1 : length(hexstring)
   if hexstring(i) ~= ' '
      bitarray(4*(i-1)+1:4*i) = bitget( hex2dec(hexstring(i)), 4:-1:1 );
   end
end
