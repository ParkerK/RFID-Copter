% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% ... a little help for ETAs.
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
% dhms = sec2dhms(seconds)
%    Translate seconds to days, hours, minutes and seconds; e.g., 86424 = 1 days, 0 hrs, 0 min,  24 sec
% 
%
%
% ***** Interface definition *****
% function dhms = sec2dhms(seconds)
%    seconds   scalar number of seconds
%    
%    dhms      SECONDS translated to [days, hours, minutes, seconds] rounded to full seconds
%
%
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
% beta 3.1   2012-08-08   arnitz      ! error in rounding operation (59.9sec was 0 min, 0 sec)
%
%
% ***** Todo *****
%
% *******************************************************************************************************

function dhms = sec2dhms(seconds)

% initialize
%     [days, hours, minutes, seconds]
dhms = zeros(1, 4);
%     round seconds to remove any rounding problems below
seconds = round(seconds);

% pick apart
%     seconds
dhms(4) = mod(seconds, 60);
seconds = fix(seconds/60);
%     minutes
dhms(3) = mod(seconds, 60);
seconds = fix(seconds/60);
%     hours, days = remainder
dhms(2) = mod(seconds, 24);
dhms(1) = fix(seconds/24);
