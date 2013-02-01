% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% sanitize settings (e.g. saturate to bounds)
%
% Warning: constraints are not checked for sanity ;-)
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
% sanitized = sanitize_setting(name, setting, [min max])
% sanitized = sanitize_setting(name, setting, constraints)
%    saturates setting to (constraints.)min <= sanitized <= (constraints.)max
% sanitized = sanitize_setting(name, setting, [min max], tol)
% sanitized = sanitize_setting(name, setting, constraints, tol)
%    saturates setting to (constraints.)min-tol <= sanitized <= (constraints.)max+tol
%
%
% ***** Function definition *****
% sanitized = sanitize_setting(name, setting, constraints)
%    name          name of setting (for warning if setting outside constraints)
%    setting       numerical value to sanitize
%    constraints   array [min,max] or struct
%       .min       enforces setting >= min
%       .max       enforces setting <= max
%    tol           (optional, default tol=0) tolerance for checks
%
%    sanitized     setting meeting constraints
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

function sanitized = sanitize_setting(name, setting, constraints, tol)


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% number of input variables
if nargin < 3
    error('Not enough input arguments.')
elseif nargin < 4
   tol = 0;
end

% if constraints is no struct, but an array => only [min, max]
if ~isstruct(constraints) && length(constraints) == 2
   constraints = struct('max', constraints(2), 'min', constraints(1));
end


% *******************************************************************************************************
% checks (allow consecutive checks)

sanitized = setting;

% out of bounds
if setting < constraints.min - tol
   if ~isempty(name); warn('Sanitizing %s (old: %g, new: %g)', name, sanitized, constraints.min); end
   sanitized = constraints.min;
end
if setting > constraints.max + tol
   if ~isempty(name); warn('Sanitizing %s (old: %g, new: %g)', name, sanitized, constraints.max); end
   sanitized = constraints.max;
end

