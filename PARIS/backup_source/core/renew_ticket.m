% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% renew system access tickets (Kerberos and AFS)
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
% [status, msg] = renew_ticket()
%    Renews Kerberos ticket and AFS authentication token, returns an error code (0: ok, 1: failed) in
%    STATUS and output / system output in MSG.
%    The function does not try to renew the ticket if the last successful renewal was within a predefined
%    timeframe. Failed renewals are repeated after a specified time (up to a maximum number of retries).
%    If the renewal failed permanently, a tolerance window is defined. A renewal will only be reported
%    "failed" if it is outside this tolerance time.
%
%
% ***** Global Variables *****
% globalsettings
%    .core.ticket_timeframe      only one ticket-renewal per timeframe (minutes)
%    .core.ticket_retries        number of immediate retries for failed renewals
%    .core.ticket_waittime       seconds waittime between immediate retries for failed renewals (seconds)
%    .core.ticket_tolerance      failed renewals within this tolerance time will not be reported (minutes)
%    .core.ticket_minrenewleft   minimum renewal time left before user interaction is requested (minutes) 
% ticket_cache     cache for renewal times, etc.
%    .lastattempt    last renewal attempt (time vector as returned by CLOCK)
%    .lastsuccess    last successful renewal (time vector as returned by CLOCK); may be empty
%    .lastfailed     last failed renewal (time vector as returned by CLOCK); may be empty
%    .numfailed      number of successively failed renewals 
%    .status         last system status report (0: ok, 1: failed)
%    .msg            last system message (output)
%
%
% ***** Interface definition *****
% function [status, msg] = renew_ticket()
%    status   1 if renewal failed, 0 otherwise
%    msg      system messages (output) if renewal failed or was skipped, empty string otherwise
%    
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
% - make system calls robust against timeouts
% = add external renew command to email output: KRB5CCNAME='echo $KRB5CCNAME' kinit -l 7d -r 30d  
% *******************************************************************************************************

function [status, msg] = renew_ticket()

% persistent ticket cache
global globalsettings;
global ticket_cache;

% does the cache exist yet?
if ~isstruct(ticket_cache)
   ticket_cache = struct('lastattempt',[], 'lastsuccess',[], 'lastfailed',[], 'numfailed',0, 'status', 0, 'msg', '');
end


% *******************************************************************************************************
% internal settings

% timestamp of maximum renewable ticket lifetime [seconds]
internalsettings.cmd_trenewable = 'date +%s -d "`klist -5 2>/dev/null | grep -A 1 krbtgt/SPSC.TUGRAZ.AT@SPSC.TUGRAZ.AT | awk ''/renew until/ {print $3" "$4}''`"';
% current timestamp [seconds]
internalsettings.cmd_tcurrent = 'date +%s';
% initialize a new ticket
internalsettings.cmd_reinit = 'kinit -l "7d" -r "30d"; aklog';

% get time (don't use tic/toc to be able to tell when a renewal failed)
currenttime = clock;


% *******************************************************************************************************
% check maximum renewable ticket lifetime

% get timestamps in seconds
[s1, msg1] = system(internalsettings.cmd_trenewable); % "renewable until"
[s2, msg2] = system(internalsettings.cmd_tcurrent); % "current time"
%     process
t1 = str2double(deblank(msg1));
t2 = str2double(deblank(msg2));

% check status / let the user enter a password if necessary
if s1~=0 || s2~=0
   err('Unable to obtain maximum ticket renewal time; system calls failed.');
else
   while t1 - t2 < globalsettings.core.ticket_minrenewleft * 60
      % inform the user that the simulation has been halted
      disp(sprintf('\n********* Ticket renewal requires user interaction. Press a button to continue. *********\n'));
      send_email('User interaction required (ticket renewal).',...
         sprintf('The maximum ticket renewal lifetime is below the threshold of %.0f minutes (%.1f days) - simulation has been halted.',...
         globalsettings.core.ticket_minrenewleft, globalsettings.core.ticket_minrenewleft/(24*60)));
      
      % user interaction
      pause; % wait for the user
      s1 = system(internalsettings.cmd_reinit); % reinitialize ticket; print messages but catch status
      
      % check again
      %     get timestamps in seconds
      [s1, msg1] = system(internalsettings.cmd_trenewable); % "renewable until"
      [s2, msg2] = system(internalsettings.cmd_tcurrent); % "current time"
      %     process
      t1 = str2double(deblank(msg1));
      t2 = str2double(deblank(msg2));
      %     output and return if successful
      if t1 - t2 >= globalsettings.core.ticket_minrenewleft * 60
         disp(sprintf('\n********* Manual ticket initialization successful *********\n'));
         ticket_cache.lastattempt = currenttime;
         ticket_cache.lastsuccess = ticket_cache.lastattempt;
         ticket_cache.lastfailed  = [];
         ticket_cache.numfailed   =  0;
         ticket_cache.status      =  0;
         ticket_cache.msg         = 'Manual ticket initialization successful.';
         status                   = ticket_cache.status;
         msg                      = ticket_cache.msg;
         return
      end
   end
end


% *******************************************************************************************************
% check if renewal is necessary

% if we had a successful renewal, this is not the first call and we are within ticket_timeframe
if ~ticket_cache.status && ~isempty(ticket_cache.lastsuccess) &&...
      etime(currenttime, ticket_cache.lastsuccess) < globalsettings.core.ticket_timeframe * 60
   msg    = 'Skipped renewal (still within skipping window)';
   status = 0;
   return
end


% *******************************************************************************************************
% renewal

% renew
for i = 1 : globalsettings.core.ticket_retries + 1
   % record this renewal
   ticket_cache.lastattempt = currenttime;
   
   % renew Kerberos ticket and AFS authentication token
   [s1, m1] = system('kinit -R');
   [s2, m2] = system('aklog');
   %     overall status
   ticket_cache.status = (s1~=0 || s2~=0);
     
   % renewal failed
   if ticket_cache.status
      warn('Kerberos ticket and/or AFS authentication token renewal failed (attempt %i of %i)', i, globalsettings.core.ticket_retries + 1);
      pause(globalsettings.core.ticket_waittime);
   else
      ticket_cache.lastsuccess = ticket_cache.lastattempt;
      break;
   end
end

% collect messages
if s1~=0 && s2~=0
   ticket_cache.msg = sprintf('%s\n******************************\n%s', deblank(m1), deblank(m2));
elseif s1~=0
   ticket_cache.msg = deblank(m1);
elseif s2~=0
   ticket_cache.msg = deblank(m2);
else
   ticket_cache.msg = '';
end

% manage tolerance window
if ticket_cache.status % failed renewal
   ticket_cache.numfailed = ticket_cache.numfailed + 1;
   if isempty(ticket_cache.lastfailed) % this was the first renewal that failed
      ticket_cache.lastfailed = currenttime;
   end
else % successful renewal: reset tolerance time
   ticket_cache.lastfailed = [];
   ticket_cache.numfailed  =  0;
end


% *******************************************************************************************************
% select return values

% if renewal failed...
if ticket_cache.status % ... and we are outside the tolerance time
   if etime(currenttime, ticket_cache.lastfailed) > globalsettings.core.ticket_tolerance * 60
      status = ticket_cache.status;
      msg    = ticket_cache.msg;
   else
      status = 0;
      msg    = sprintf('Within tolerance time for failed renewals (%i failed)', ticket_cache.numfailed);
   end

% everything ok
else
      status = ticket_cache.status;
      msg    = ticket_cache.msg;
end
