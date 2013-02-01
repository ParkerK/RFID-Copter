% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% bit array to settings
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
% opt = bit2opt(bitarray)
%    Will return the unsigned integer representation of BITARRAY (MSB = bitarray(1)). 
%    Idential to bit2opt(bitarray, {0,1,2,...,2^(length(bitarray)-1)}).
% opt = bit2opt(bitarray, options)
%    The unsigned integer representation of BITARRAY is used as index (+1) in options.
%    Returns options(integer(bitarray)+1).
%
%
% ***** Interface definition *****
% setting = bit2opt(bitarray, options)
%    bitarray   binary array to translate to decimal or option. Interpreted as unsigned integer.
%               (e.g. [1,1,1,0]=14)
%    options    (optional) cell array containing options to choose from
%               If provided, bit2set will use bitarray as index in options and and return the cell
%               bitarray is pointing at. options{1} corresponds to bitarray=[0,...,0]
%
%    opt        integer(bitarray) or options(integer(bitarray)+1).
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

function opt = bit2opt(bitarray, options)

% *******************************************************************************************************
% input parameter checks

% number of input parameters
if nargin < 1
   error('Not enough input arguments.');
elseif nargin == 1
   options = {};
end

% empty bitarray
if isempty(bitarray)
   opt = 0;
   return
end

% options has to be cell array
if ~iscell(options)
   error('options has to be a cell array');
end

% make bitarray column vector
bitarray = bitarray(:);


% *******************************************************************************************************
% "the function"

% bit array => number
opt = sum(bitarray .* 2.^[length(bitarray)-1:-1:0]');

% if options are provided: number -> option
if ~isempty(options)
   % index outside options?
   if opt+1 > length(options)
      warn('Binary index exceeds length of options. Returning index number instead.');
   else
      opt = options{opt+1};
   end
end
