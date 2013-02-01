% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% Everything needed for ETA estimates in loops assuming constant execution time for all passes.
%
% Attention: 
% .) Matlab's stopwatch timer has to be initialized using TIC before calling this function for the first
%    time.
% .) Initialize before starting the loop and call this function at the **beginning** of a new loop pass.
%       eta_forloop = get_eta();
%       for 1 = 1 : 100
%          eta_for = get_eta(eta_forloop, i, 100);
%          ...
%       end
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
% eta123 = get_eta()
%    Initialize a new eta struct with the default values listed below.
% eta123 = get_eta(eta123, loop_now, loop_len)
%    Update the eta struct ETA123 using the current loop counter LOOP_NOW and the total number of loop
%    passes LOOP_LEN.
%    
%
%
% ***** Interface definition *****
% function eta = get_eta(eta, loop_now, loop_len)
%    eta              struct containing timing information; default values are listed in parenthesis
%       .ticID           identifier of the TIC command for this ETA structure (allows nesting of ETAs)
%       .sec_past        seconds past since the creation of this struct (0)
%       .sec_left        seconds left until end of loop (NaN)
%       .time_past       [days, hours, minutes, seconds] past since the creation of this struct (NaNs)
%       .time_left       [days, hours, minutes, seconds] left until end of loop (NaNs)
%       .datenum_start   time of generation for this struct as a serial date number (NOW)
%       .datenum_stop    ETA as serial date number (NaN)
%    loop_now         current loop index
%    loop_len         maximum loop index (number of loop passes)
%    
%    eta    newly created eta struct or updated struct based on LOOP_NOW and LOOP_LEN
%
%
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
% beta 3.1   2012-06-06   arnitz      ! replaced sec_start by tic ID number; allows nesting of eta cmds
%
%
% ***** Todo *****
%
% *******************************************************************************************************

function eta = get_eta(eta, loop_now, loop_len)

% called without parameters: generate a new eta struct
if nargin == 0
   eta = struct('ticID', tic, 'sec_past',0, 'sec_left',NaN,...
      'time_past',zeros(1,4), 'time_left',nan(1,4), 'datenum_start',now, 'datenum_stop',NaN);

% otherwise: calculate ETA (if possible)
else
   % get times [seconds]
   eta.sec_past = toc(eta.ticID); % up to now
   eta.sec_left = eta.sec_past / (loop_now-1) * (loop_len - loop_now + 1);
   % get times [days, hours, minutes, seconds]
   eta.time_past = sec2dhms(eta.sec_past);
   eta.time_left = sec2dhms(eta.sec_left);
   % get ETA
   eta.datenum_stop = now + eta.sec_left/86400;
end


