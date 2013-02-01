% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% get stack level (counting only primary functions)
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
% level = get_stacklevel()
%     Returns the number of primary functions in stack including the calling function
%     (not including get_stacklevel). Does not count core functions as primary function (they are
%     supposed to be in the background). Uses the function ISCOREFUNCTION do determine if a function 
%     is a core function.
% level = get_stacklevel(n)
%     Like get_stacklevel(), but skipping the first n entries. Use get_stacklevel(-1) to include this
%     function as well.
%
%
% ***** Interface definition *****
% level = get_stacklevel(n)
%    n       omit first n entries in stack
%
%    level   number of primary functions in stack
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

function level = get_stacklevel(n)

% default setting for skip
if nargin == 0
   n = 0;
end

% check if n is numeric
if ~isnumeric(n)
   error('Parameter n (number of stackentries to skip) has to be numeric.');
end
  
% get stack
stack = dbstack;
stack(1+n) = []; % remove this + n functions (dbstack(1) does not really work with Matlab < 2007a)

% determine level, ignore nonprimary functions and core functions
level = 0;
for i = 2 : length(stack)
   if strcmp(stack(i).file(1:end-2), stack(i).name) &&... % if this is a primary function or script
         ~iscorefunction(stack(i).name)                   % this is not a core function
      level = level + 1;
   end
end

