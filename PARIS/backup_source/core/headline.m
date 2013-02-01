% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% logging system - print a headline message
% Used for: headlines in main modules (can be switched off)
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
% headline(message, a1, a2, ...)
%    Print MESSAGE in the command window. Can handle sprintf-like arguments directly.
%
% Can be activated/deactivated via globalsettings.logging.headlines (default: active)
%
%
% ***** Global Variables *****
% globalsettings
%    .logging.headlines   print headlines (true/false)?
%
%
% ***** Interface definition *****
% function headline(text, a1, a2, ...)
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
%
% *******************************************************************************************************

function headline(varargin)

% global settings
global globalsettings;

% check if warnings should be printed
if isstruct(globalsettings) && isstruct(globalsettings.logging) &&...
      isfield(globalsettings.logging, 'headlines') && ~globalsettings.logging.headlines     
   return
end

% first cell of varargin has to be string (sprintf string)
varargin{1} = sprintf('%s', varargin{1});

% issue the warning 
disp(sprintf(varargin{:}));
