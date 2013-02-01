% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - finds and returns the indices of zero crossings (quick, dirty and nonperformant)
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
% version = findzeros()
%    Just returns the version number (string).
% zerocrossings = findzeros(signal)
%    Returns all indices next to zero crossings in SIGNAL. If SIGNAL has no zero crossings or is a
%    scalar, an empty vector is returned. Matrices will be serialized (linear indexing).
%
%
% ***** Function definition *****
% function zerocrossings = findzeros(signal)
%    signal          signal vector (or matrix)
%
%    zerocrossings   all indices of zero crossings in signal
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
% = more performant implementation (massive speed gain for self test)
%   
% *******************************************************************************************************

function zerocrossings = findzeros(signal)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   zerocrossings = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% number of input variables
if nargin < 1
    criterr('Not enough input arguments.')
end

% length of input vector
if isempty(signal)
    zerocrossings = [];
    return
end

% if signal is not a vector => serialize
if ~isvector(signal)
   signal = signal(:);
end


% *******************************************************************************************************
% quick and dirty zero crossing detection

zerocrossings = [];
for i = 1 : length(signal) - 1
   % exact zero
   if signal(i) == 0
      zerocrossings = [zerocrossings; i];
   % falling edge || rising edge
   elseif (signal(i) > 0 && signal(i+1) < 0) || (signal(i) < 0 && signal(i+1) > 0)
      if abs(signal(i)) < abs(signal(i+1))
         zerocrossings = [zerocrossings; i];
      else
         zerocrossings = [zerocrossings; i+1];
      end
   end
end   
