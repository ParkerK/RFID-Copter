% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% compare two floats for "approximate" equality (array/scalar <-> array/scalar)
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
% eq = approxequal(a, b, thr)
%    Returns true if A and B are approximately equal, i.e., if the rel. difference is smaller than thr.
% eq = approxequal(a, b)
%    Uses a default threshold of 1000*eps.
% [eq, err] = approxequal(a, b, thr)
%    Also returns the relative difference (error) between A and B.
%
%
% ***** Interface definition *****
% function eq = approxequal(a, b, thr)
%    a     first value to compare
%    b     second value to compare
%    thr   threshold for relative difference
%
%    eq    true if A and B are approximately equal, false if not
%    err   
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

function [eq, err] = approxequal(a, b, thr)

% default threshold?
if ~exist('thr', 'var')
   thr = 1000*eps;
end

% in case of arrays
a = a(:);
b = b(:);

% length check
if length(a) ~= 1 && length(b) ~= 1 && length(a) ~= length(b)
   err = NaN;
   eq  = false;
end

% relative difference
if all(a == 0) && all(b == 0)
   err = 0;
elseif max(abs(a)) ~= 0
   err = abs(a-b) ./ max(abs(a));
elseif max(abs(b)) ~= 0
   err = abs(a-b) ./ max(abs(b));
else
   err = Inf; % this should not happen
end

% approximately equal?
eq = err < thr;
   

