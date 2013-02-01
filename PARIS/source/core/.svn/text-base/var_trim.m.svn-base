% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% variance of a vector ignoring leading and trailing zeros
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
% Calculates the variance of a vector (arrays are transformed to vectors) ignoring leading and trailing
% zeros. Otherwise the behavior is identical to that of VAR.
%
%
% ***** Interface definition *****
% see Matlab help for VAR
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

function v = var_trim(x, flag, dim)

% a quick check
if nargin < 1
   err('Not enough input parameters.');
elseif nargin > 3
   err('Too many input parameters.');
end

% flatten arrays
x = x(:);

% get indices of first and last nonzero element
ind1 = find(x~=0, 1, 'first');
ind2 = find(x~=0, 1, 'last');

% calculate variance
if nargin == 1
   v = var(x(ind1:ind2));
elseif nargin == 2
   v = var(x(ind1:ind2), flag);
else
   v = var(x(ind1:ind2), flag, dim);
end
