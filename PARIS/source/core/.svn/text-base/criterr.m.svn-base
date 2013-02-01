% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% throw critical error
% Used for: program errors (e.g. function in undefined state because previous checks failed)
%
% ATTENTION:
% Do not call criterr from a function that does not support version calls, because it tries to get the
% version of the calling function. If such a call is not supported, this call may result in
% malfunction or data loss.
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
% criterr()
%    Like criterr('Unspecified error')
% criterr(message, a1, a2, ...)
%    Throws an error with identifier "SPR_sim:[throwing primary function]:critical" and additional text
%    "Critical error". Version information for all functions involved is added to stack.version (see
%    DBSTACK for further details on stack). Can handle sprintf-like arguments directly (like ERROR does).
%    An add. stack output can be (de)activated via globalsettings.logging.exceptions (default: active)
%
%
% ***** Global Variables *****
% globalsettings
%    .logging.exceptions   print (additional) error messages (true/false)?
%    .core.debug           switch to debug mode in case of an error
%
%
% ***** Interface definition *****
% function criterr(message, a1, a2, ...)
%    see Matlab help for "error('message')" or "error('message', a1, a2, ...)" for details
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
% - Matlab displays "error using ==> err at" when called from command line (fix necessary ?)
%
% *******************************************************************************************************

function criterr(varargin)

% global variables
global globalsettings;

% if varargin is empty (no message has been provided)
if isempty(varargin)
   varargin = {'Unspecified error'};
end

% browse through dbstack and add version numbers (do not include this function)
stack = dbstack;
stack(1) = []; % workaround for dbstack(1) for Matlab R2006a,b
for i = 1 : length(stack)
   status = mfile_status(stack(i).file); % get status of this file
   stack(i).version = status.version;
end

% check for varargin{1}: not a string
% and varargin{:}: not a valid sprintf argument
% ... borrowed identifier from sprintf(1)
% ... unfortunately Matlab displays "error using ==> err at" when called from command line
if ~ischar(varargin{1}) || isempty(sprintf(varargin{:}))
   error(struct('message', 'Invalid format.', 'identifier', 'MATLAB:badformat_mx', 'stack', stack));
end
   
% first cell of varargin is a string => we can easily add something
varargin{1} = sprintf('Critical error: %s', varargin{1});

% create message struct
if isempty(stack) % called from command line (unlikely, but still it shouldn't produce an error)
   error('Cannot be called from command line.');
else
   message_struct = struct('message', sprintf(varargin{:}), ...
      'identifier', sprintf('SPR_sim:%s:critical', stack(1).file(1:end-2)), ...
      'stack', stack);
end

% create error message (for output and email)
error_msg = sprintf('%s\n\nIdentifier:\n%s', message_struct.message, message_struct.identifier);
if ~isempty(length(message_struct.stack))
   error_msg = sprintf('%s\n\nStack Trace:', error_msg);
   for i = length(message_struct.stack) : -1 : 1
      error_msg = sprintf('%s\n%s-> %s %s, line %3.0f (%s)', error_msg, blanks(length(message_struct.stack)-i), ...
         message_struct.stack(i).file, message_struct.stack(i).version,...
         message_struct.stack(i).line, message_struct.stack(i).name);
   end
end

% debug mode (if requested)
if isstruct(globalsettings) && isstruct(globalsettings.core) &&...
      isfield(globalsettings.core, 'debug') && globalsettings.core.debug
   dbstop if error
   error_msg = sprintf('%s\n\nEntered debug mode.', error_msg);
end

% display message (if requested)
if isstruct(globalsettings) && isstruct(globalsettings.logging) &&...
      isfield(globalsettings.logging, 'exceptions') && globalsettings.logging.exceptions
   disp(error_msg);
%    disp(sprintf('\n\n'));
%    diary off; % leave inside if structure (otherwise diary may be stopped during selftests)
end

% send an email
% send_email('Simulation terminated abnormally with critical error', error_msg);

% issue the error
error(message_struct)
