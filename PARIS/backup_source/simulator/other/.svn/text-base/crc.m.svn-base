% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - general CRC implementation
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
% version = crc()
%    Just returns the version number (string).
% crc_value = crc(data, 5)
%    Performs the CRC-5 encoding according to EPCglobal class-1 gen-2.
% crc_value = crc(data, 16)
%    Performs the CRC-16 encoding according to EPCglobal class-1 gen-2.
% crc_value = crc(data, poly)
%    Performs standard CRC encoding with SREG initialization=0x0 and returns the non-inverted remainder 
% crc_value = crc(data, poly, init, invert)
%    Performs the CRC encoding with SREG initialization=init and returns the non-inverted remainder, if
%    invert=false and inverted remainder if invert=true. 
%
%
% ***** Interface definition *****
% function crc_value = crc(data, poly, init, invert)
%    data     input data to encode a_m x^m + a_(m-1) x^(m-1) + ... + a_0 => [a_m, ... , a_0]
%    poly     CRC polynomial       b_m x^m + b_(m-1) x^(m-1) + ... + b_0 => [b_k, ... , b_0]
%             example: polynomial x^5+x^3+1 => poly=[1,0,1,0,0,1]
%             or type of CRC (only EPCglobal class-1 gen-2 CRC-5 and CRC-16 supported)
%    init     (optional) SREG initialization (preload SREG with init)
%             Assumes that bits are shifted in from the left (standard signal flow), thus
%             init(1) corresponds to the rigthmost bit and thus to data(1).
%    invert   (optional) returns the inverted remainder if true, non-inverted if false
%
%    crc      CRC result (remainder)
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

function crc_value = crc(data, poly, init, invert)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   crc_value = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% CRC-5 generator according to EPCglobal class-1 gen-2 p. 86
internalsettings.poly5    = [1;0;1;0;0;1]; % x^5+x^3+1
internalsettings.init5    = [0;1;0;0;1];   % preload SREGs with C[4:0]=01001b and
internalsettings.invert5  = 0;             % do not invert

% CRC-16 generator according to EPCglobal class-1 gen-2 p. 86
internalsettings.poly16   = [1;0;0;0;1;0;0;0;0;0;0;1;0;0;0;0;1]; % x^16+x^12+x^5+1
internalsettings.init16   = ones(16,1);                          % preload SREGs with 0xFFFF and
internalsettings.invert16 = 1;                                   % invert


% *******************************************************************************************************
% input parameter checks

% number of input parameters
if nargin == 2
   if length(poly) == 1 % EPCglobal CRC-5 or CRC-16
      switch poly
         case 5
            poly   = internalsettings.poly5;
            init   = internalsettings.init5;
            invert = internalsettings.invert5;
         case 16
            poly   = internalsettings.poly16;
            init   = internalsettings.init16;
            invert = internalsettings.invert16;
         otherwise
            err('Only CRC-5 or CRC-16 allowed for CRC(data, type). See help crc for details.');
      end
   elseif length(poly) > 1 % standard CRC with init=zeros and noninverted remainder
      init   = zeros(length(poly)-1,1);
      invert = 0;
   else
      err('Length of polynomial vector is zero.')
   end
elseif nargin == 3
   criterr('Not enough input arguments.');
elseif nargin < 2
   criterr('Not enough input arguments.');
end

% empty vectors
if isempty(data)
   err('Length of data vector is zero.')
end

% length of init vector has to be length of polynomial
if length(init) ~= length(poly)-1
   err('init has to have same length as poly-1.')
end

% nonbinary alphabet
if min(data) < 0 || max(data) > 1
   err('data has non-binary alphabet.')
end
if min(poly) < 0 || max(poly) > 1
   err('poly has non-binary alphabet.')
end
if min(init) < 0 || max(init) > 1
   err('init has non-binary alphabet.')
end

% invalid polynomial
if poly(1) ~= 1
   err('poly(1), i.e. x^k must not be zero.')
end

% create column vectors
data = data(:);
poly = poly(:);


% *******************************************************************************************************
% CRC implementation (with SREG init 0xFF...)

% length of polynomial
k = length(poly) - 1; % (x^k + x^(k-1) + ... + 1)

% prepare data
data      = [data; zeros(k,1)];     % append k bits ("x^k for data")
data(1:k) = xor( data(1:k), init ); % preload with init vector is like selective inversion

% CRC calculation (straightforward polynomial division)
for i = 1 : length(data)-k
   if data(i) == 1
      data(i:i+k) = xor(data(i:i+k), poly);
   end
   
end

% truncate to length of CRC and invert
if invert
   crc_value = xor( data(end-k+1:end), ones(k,1) );
else
   crc_value = data(end-k+1:end);
end


