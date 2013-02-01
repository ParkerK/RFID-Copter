% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% channel - (antenna) directivity
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
% version = channel_directivity()
%    Just returns the version number (string).
% gain = channel_directivity(settings)
%    Returns the combined gain factor of TX and RX antenna for one channel according to SETTINGS.
% [gain, dir_tx, dir_rx] = channel_directivity(settings)
%    Returns the selected angles (azimuth and elevation) in both RX and TX gain patterns along with the
%    combined gain (useful mostly for debugging).
%
%
% ***** Interface definition *****
% function [gain, dir_tx, dir_rx] = channel_directivity(settings)
%    settings   struct containing all necessary parameters
%       .txant    filename of TX antenna characteristic (isotropic if empty)
%       .txrot    rotation of TX antenna [azimuth, elevation] in degree ([0,0]: maximum in positive x-axis)
%       .rxant    filename of RX antenna characteristic (isotropic if empty)
%       .rxrot    rotation of RX antenna [azimuth, elevation] in degree ([0,0]: maximum in positive x-axis)
%       .dir_tx   direction of channel [azimuth, elevation] in degree (from TX antenna)
%       .dir_rx   (optional; if the path RX->TX is NLOS) like .dir_tx, from a fictive TX antenna
%                 resulting in a direct path (will be used to calculate RX antenna gain)
%
%   gain        combined gain factor of RX and TX antenna
%   dir_tx      angles [azimuth, elevation] in degree in the TX gain pattern
%   dir_rx      angles [azimuth, elevation] in degree in the RX gain pattern
%
%
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
%      3.1   2012-05-15   arnitz      ~ also returns angles in RX/TX patterns (nice for debugging) 
%
%
% ***** Todo *****
% - third rotation dimension: [azimuth, elevation, rotation]
%
% *******************************************************************************************************


function [gain, dir_tx, dir_rx] = channel_directivity(settings)
version = 'beta 3.1';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   gain = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks / prepare parameters

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'txant', 'txrot', 'rxant', 'rxrot', 'dir_tx'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% optional parameter settings.dir_rx
if ~isfield(settings, 'dir_rx')
   settings.dir_rx = settings.dir_tx;
end

% check if all angles contain azimuth and elevation
if (size(settings.txrot, 2) ~= 2) || (size(settings.rxrot, 2) ~= 2) ||...
      (size(settings.dir_tx, 2) ~= 2) || (size(settings.dir_rx, 2) ~= 2)
   err('All angles in settings have to contain azimuth and elevation ([az, el]).');
end

% check size of settings.dir_rx and settings.dir_tx
if (numel(settings.dir_tx) ~= numel(settings.dir_rx)) ||...
      (ndims(settings.dir_tx) ~= ndims(settings.dir_rx)) || any( size(settings.dir_tx) ~= size(settings.dir_rx) )
   err('Size of settings.dir_tx and settings.dir_rx have to match.');
end


% *******************************************************************************************************
% calculate TX/RX antenna gain and combine

% "rotate antennas" [azimuth, elevation]; dir_rx, dir_tx: angles in antenna gain diagrams
%     ... directions are measured from transmitter => receiver rotated to match
settings.dir_tx = mod(settings.dir_tx - repmat(settings.txrot,            size(settings.dir_tx,1), 1), 360);
settings.dir_rx = mod(settings.dir_rx - repmat(settings.rxrot + [0, 180], size(settings.dir_rx,1), 1), 360); % +[0,180] = "on the other side"
%     adapt to standard sphere coordinates (azimuth 0..360, elevation 0..180 degree)
settings.dir_tx = mod(settings.dir_tx, 360);
settings.dir_tx(settings.dir_tx(:,2)>180,1) =   mod(settings.dir_tx(settings.dir_tx(:,2)>180,1) - 180, 360);
settings.dir_tx(settings.dir_tx(:,2)>180,2) = 360 - settings.dir_tx(settings.dir_tx(:,2)>180,2);
settings.dir_rx = mod(settings.dir_rx, 360);
settings.dir_rx(settings.dir_rx(:,2)>180,1) =   mod(settings.dir_rx(settings.dir_rx(:,2)>180,1) - 180, 360);
settings.dir_rx(settings.dir_rx(:,2)>180,2) = 360 - settings.dir_rx(settings.dir_rx(:,2)>180,2);

% transmitter antenna
if isempty(settings.txant)
   gain_tx = ones(size(settings.dir_tx,1), 1);
else
   % load data
   txdata = loadmat(settings.txant);
   % interpolate: 2x2D or 3D
   if txdata.full3d
      [az, el] = meshgrid(txdata.az, txdata.el);
      gain_tx  = interp2(az, el, txdata.gain', settings.dir_tx(:,1), settings.dir_tx(:,2), 'linear', NaN);
      if any(isnan(gain_tx))
         err('TX antenna gain contains NaNs. This is likely a numerical issue (angle slighty outside 0..360).');
      end
   else
      gain_tx = ...
         interp1(txdata.az, txdata.az_gain, settings.dir_tx(:,1), 'linear', 'extrap') .* ...
         interp1(txdata.el, txdata.el_gain, settings.dir_tx(:,2), 'linear', 'extrap'); % do not use spline: discontinuities around gain==0
   end
end

% receiver antenna
if isempty(settings.rxant)
   gain_rx = ones(size(settings.dir_rx,1), 1);
else
   % load data
   rxdata = loadmat(settings.rxant);
   % interpolate: 2x2D or 3D
   if rxdata.full3d
      [az, el] = meshgrid(rxdata.az, rxdata.el);
      gain_rx  = interp2(az, el, rxdata.gain', settings.dir_rx(:,1), settings.dir_rx(:,2), 'linear');
      if any(isnan(gain_tx))
         err('RX antenna gain contains NaNs. This is likely a numerical issue (angle slighty outside 0..360).');
      end
   else
      gain_rx = ...
         interp1(rxdata.az, rxdata.az_gain, settings.dir_rx(:,1), 'linear', 'extrap') .* ...
         interp1(rxdata.el, rxdata.el_gain, settings.dir_rx(:,2), 'linear', 'extrap'); % do not use spline: discontinuities around gain==0
   end
end

% combine gains
gain = gain_tx .* gain_rx;

% als return angles; might be of interest in debugging
dir_tx = settings.dir_tx;
dir_rx = settings.dir_rx;

