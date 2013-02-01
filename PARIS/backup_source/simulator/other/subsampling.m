% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - track frequency components in case of subsampling
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
% version = subsampling()
%    Just returns the version number (string).
% fi_rs = subsampling(fi, frxs)
%    Returns the location of the frequency components FI after sampling with receiver sampling frequency
%    FRXS.
%   
%
%
% ***** Interface definition *****
% function fi_rs = subsampling(fi, frxs)
%    fi      location of frequency components before sampling
%    frxs    receiver sampling frequency in Hz
%
%    fi_rs   lovation of frequency components after sampling
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


function fi_rs = subsampling(fi, frxs)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   fi_rs = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% track frequency components

fi_rs = fi - frxs*floor(fi/frxs);
fi_rs(fi_rs > frxs/2) = fi_rs(fi_rs > frxs/2) - frxs;

