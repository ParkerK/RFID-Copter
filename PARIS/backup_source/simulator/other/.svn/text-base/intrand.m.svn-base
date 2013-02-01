% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - generate uniformly distributed integer random variables
%
% Note that "integer" here means the set of numbers (natural numbers including their negatives and zero)
% and not the Matlab type.
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
% version = intrand()
%    Just returns the version number (string).
% out = intrand(size, bounds)
% out = intrand(bounds)
%    Returns an array of size N with discrete integer random variables uniformly distributed between 
%    BOUNDS(1) and BOUNDS(2) using Matlab's builtin function RAND. Like for RAND, the size parameter can
%    have different formats (m,n,p,... or [m,n,p,...] to create an array of size m-by-n-by-p-by...). 
%    For example (uniform distribution between a and b):
%       intrand(2,3,1,[a,b])  creates a 2-by-3-by-1 array; identical to intrand([2,3,1],[a,b])
%       intrand(3,[a,b])      creates a 3-by-3 array; identical to intrand(3,3,[a,b])
%       intrand([a,b])        creates a scalar; identical to intrand(1,[a,b])   
%
%
% ***** Interface definition *****
% function out = intrand(size, bounds)
%    size     desired size of array OUT (see Matlab help: RAND)
%    bounds   bounds of uniform distribution [lower, upper] (have to be integers ...-1,0,1,2,...)
%
%    out      array with integer samples uniformly distributed between BOUNDS(1) and BOUNDS(2)
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
% ? switch to randi (available since Matlab 2008b)
%
% *******************************************************************************************************

function out = intrand(varargin)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   out = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% check/prepare parameters and generate random vector/matrix

% "demultiplex" input
if nargin == 1
   n = 1;
else
   n = cell2mat(varargin(1:end-1));
end
b = varargin{end};

% check if bounds are integer (in natural numbers)
if ~all(round(b) == b)
   err('Bounds for integer random variable have to be integer.');
end

% create vector
% ... eps-0.5 and 1-2*eps is necessary to make boundaries equally probable 
%     ( bounds of b(1)*rand(b(2)-b(1)) have p/2, where p is the probability within the bounds )
out = b(1) + round( eps-0.5 + rand(n)*(b(2)-b(1)+1-2*eps) );
