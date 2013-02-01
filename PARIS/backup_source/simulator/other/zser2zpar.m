% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - parallel <-> serial equivalent impedance conversion (for one w)
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
% version = zser2zpar()
%    Just returns the version number (string).
% zout = zser2zpar(zin, 's2p')
%    Converts serial equivalent impedance ZIN to parallel equivalent impedance ZOUT.
% zout = zser2zpar(zin, 'p2s')
%    Converts parallel equivalent impedance ZIN to serial equivalent impedance ZOUT.
%
%
% ***** Interface definition *****
% function zout = zser2zpar(zin, mode)
%    zin    complex input impedance  (Zser or Zpar depending on MODE)
%    mode   serial->parallel (s2p) or parallel->serial (p2s) conversion {'s2p', 'p2s'} 
%
%    zout   complex output impedance (Zser or Zpar depending on MODE)
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

function zout = zser2zpar(zin, mode)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   zout = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% conversion

switch lower(mode)
   case 'p2s'
      rp = real(zin);
      xp = imag(zin);
      zout = complex(rp.*xp.^2./(rp.^2+xp.^2), rp.^2.*xp./(rp.^2+xp.^2));
   case 's2p'
      rs = real(zin);
      xs = imag(zin);
      rp = rs + xs.^2./rs;
      xp = rp.*sqrt(rs)./sqrt(rp-rs) .* sign(xs);
      zout = complex(rp, xp);
   otherwise
      err('Unsupported mode "%s". Only "p2s" (parallel -> serial) and "s2p" (serial -> parallel) possible.', lower(mode));
end


