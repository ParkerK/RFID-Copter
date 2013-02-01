% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% logging system - issue a message
% Used for: every message that is not a warning (e.g. version information of characteristics)
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
% msg(message, a1, a2, ...)
%    Print a text in the command window Can handle sprintf-like arguments directly.
%
% Can be activated/deactivated via globalsettings.logging.messages (default: active)
%
%
% ***** Global Variables *****
% globalsettings
%    .logging.messages   print messages (true/false)?
%
%
% ***** Interface definition *****
% function msg(message, a1, a2, ...)
%    see Matlab help for sprintf for details
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
% - verbosity level to be able to display certain messages without displaying all of them
%
% *******************************************************************************************************

function msg(varargin)

% global settings
global globalsettings;

% check if warnings should be printed
if isstruct(globalsettings) && isstruct(globalsettings.logging) &&...
      isfield(globalsettings.logging, 'messages') && ~globalsettings.logging.messages     
   return
end

% first cell of varargin has to be string ... add indentation and leading text
varargin{1} = sprintf('%s%s', blanks(3*(get_stacklevel()+1)), varargin{1});

% issue the warning 
disp(sprintf(varargin{:}));
