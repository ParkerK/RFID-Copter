% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% logging system - issue critical warning
% Used for: settings that might cause an error or data loss/corruption
%           (e.g. manual override of housekeeping functions)
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
% critwarn(message, a1, a2, ...)
%    Issues a warning with identifier "SPR_sim:[throwing primary function]:critical" and additional text
%    "WARNING". Can handle sprintf-like arguments directly (like WARNING does).
%
% Can be activated/deactivated via globalsettings.logging.warnings (default: active)
%
%
% ***** Global Variables *****
% globalsettings
%    .logging.warnings   print warning messages (true/false)?
%
%
% ***** Interface definition *****
% function critwarn(message, a1, a2, ...)
%    see Matlab help for "warning('message')" or "warning('message', a1, a2, ...)" for details
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

function critwarn(varargin)

% global settings
global globalsettings;

% check if warnings should be printed
if isstruct(globalsettings) && isstruct(globalsettings.logging) &&...
      isfield(globalsettings.logging, 'warnings') && ~globalsettings.logging.warnings     
   return
end

% first cell of varargin has to be string ... add indentation and leading text
varargin{1} = sprintf('%s=> WARNING: %s', blanks(3*(get_stacklevel()+1)), varargin{1});

% issue the warning 
disp(sprintf(varargin{:}));
