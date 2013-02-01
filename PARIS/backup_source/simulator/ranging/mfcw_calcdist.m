% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% ranging - MFCW distance estimate the (noncomplex/complex) MFCW derivation
%
% Note: Assumes \omega_0 = 0, hence SETTINGS.FI(i) = \omega_i, but c_mi(1) corresponds to \omega_0!
% Note: 2*n-FCW does not reuse carriers (iid. theory applies), (2*n-1)-FCW uses all carrier => overdef!
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
% version = mfcw_calcdist()
%    Just returns the version number (string).
% dist = mfcw_calcdist(c_mi, c_im, settings)
%    Returns distance estimates out of the averaged component amplitudes C_MI (lower sideband:
%    fi-modulation) and C_IM (upper sideband fi+modulation). The results of i-frequency CW ranging is
%    stored in DIST{i}. Systematic errors in C_MI and C_IM have to be removed beforehand.
%   
%
%
% ***** Interface definition *****
% function dist = mfcw_calcdist(c_mi, c_im, settings)
%    c_mi   complex amplitudes of lower sidebands ([0, fi]-modulation)
%    c_im   complex amplitudes of upper sidebands ([0, fi]+modulation)
%    settings   struct containing settings
%       .nc        number of secondary carriers to consider (.nc has to be <= length(.fi)!)
%       .fi        secondary secondary carrier offset frequencies in Hz
%       .c_ord     order in which to combine the carriers (USB and LSB); leave empty for default
%                  permutation {[], [1,2], [1,2; 1,3], [1,2; 3,4], [1,2; 3,4; 1,5], [1,2; 3,4; 5,6]}
%       .c         speed of light in m/s
%
%    dist   one distance estimate for each permutation in SETTINGS.C_ORD (note that dist{1}=NaN)
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
% - detect and remove outliers
%
% *******************************************************************************************************


function dist = mfcw_calcdist(c_mi, c_im, settings)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   dist = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% determine order of phase differences (carrier indices)
internalsettings.c_ord{2} = [1,2];
internalsettings.c_ord{3} = [1,2; 1,3];
internalsettings.c_ord{4} = [1,2; 3,4];
internalsettings.c_ord{5} = [1,2; 3,4; 1,5];
internalsettings.c_ord{6} = [1,2; 3,4; 5,6];


% *******************************************************************************************************
% input parameter checks / prepare parameters

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'nc', 'fi', 'c_ord', 'c'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% no cw-ranging information whatsoever...
if settings.nc == 0
   msg('No MFCW secondary carriers. Is this a mistake?');
   dist = {};
   return
end

% empty carrier order
if isempty(settings.c_ord)
   settings.c_ord = internalsettings.c_ord;
end


% *******************************************************************************************************
% estimate distance

% estimate phase shifts (all sidebands)
for i = 2 : settings.nc + 1 % 2, 3, 4, ..., frequency CW
   dphi{i} = phase_diff(c_mi, c_im, settings.c_ord{i});
end

% extract frequencies from settings and add f0 = 0 (mfcw_addseccarriers assumes f0=0)
wi = 2*pi*[0; settings.fi(:)];

% phase shift -> distance
for i = 2 : settings.nc + 1 % 2, 3, 4, ..., frequency CW
   % individual distance estimates
   domega = repmat( wi(settings.c_ord{i}(:,2))' - wi(settings.c_ord{i}(:,1))', 1, 2);
   dist_hat = dphi{i} * settings.c ./ (2*domega); % LSB, USB symmetrical
   % try to remove ambiguities
   lambda = 2*pi*settings.c ./ abs(domega);
   dist_hat(dist_hat < -lambda/4) = dist_hat(dist_hat < -lambda/4) + lambda(dist_hat < -lambda/4)/2;
   dist_hat(dist_hat < -lambda/8) = dist_hat(dist_hat < -lambda/8) + lambda(dist_hat < -lambda/8)/4;
   % average
   dist{i} = mean(dist_hat);
end

end



% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% phase differences between USB (LSB) of carriers
function D = phase_diff(c_mi, c_im, c_order)
for i = 1 : size(c_order, 1)
   D(i                ) = angle( c_mi(c_order(i,1)) ./ c_mi(c_order(i,2)) );
   D(i+size(c_order,1)) = angle( c_im(c_order(i,1)) ./ c_im(c_order(i,2)) );
end
end
