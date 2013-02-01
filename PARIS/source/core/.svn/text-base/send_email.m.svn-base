% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% send an email (if available)
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
% send_email(subject, body)
%    Sends an email to the address specified in GLOBALSETTINGS if an automated email system is 
%    available (currently only supported on properly configured Unix systems). Note that quotation marks
%    (") in SUBJECT and BODY will be replaced by ticks (');
%
%
% ***** Global Variables *****
% globalsettings
%    .core.mailto        email address (no mail sent if empty)
%    .core.subj_prefix   prefix for subject line to identify the mail in email clients
%    .misc.hostname      name of host we're running on
%    .misc.pid           own process ID
%
%
% ***** Interface definition *****
% function send_email(subject, body)
%    subject   message subject (string)
%    body      message body (string)
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
% - add a field "which simulation" to body (if possible)
% ? switch to sendmail
%
% *******************************************************************************************************

function send_email(subject, body)

% global variables
global globalsettings;

% check input
if isempty(subject) || ~ischar(subject)
   error('Subject is empty or not a string.');
end
%     subject and body must not contain quotation marks (shell command below)
if any(subject == '"')
   subject(subject == '"') = '''';
end
if any(body == '"')
   body(body == '"') = '''';
end

% send the email
if isunix && ~isempty(globalsettings.core.mailto)
   system(sprintf('echo "%s" | mail --subject="%son %s (pid %i): %s" %s', body, globalsettings.core.subj_prefix,...
      globalsettings.misc.hostname, globalsettings.misc.pid, subject, globalsettings.core.mailto));
end
