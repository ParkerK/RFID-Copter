% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% channel - large-scale (wideband, time-invariant)
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
% version = channel_large()
%    Just returns the version number (string).
% [delay_s, gain] = channel_large(settings)
%    Returns delay and gain factor for the large-scale channel between receivers and transmitters
%    according to SETTINGS. Currently only a log-distance model is supported. Note that distances are 
%    not checked for (geometric) sanity! Supports up to two-dimensional distance settings 
%    (all distances transmitter -> receiver can be processed at once).
%    
%
%
% ***** Interface definition *****
% function [delay_s, gain] = channel_large(settings)
%    settings   struct containing all necessary parameters
%       .type      type of large-scale model {'log-dist'}
%       .dist      distance/distances in m (supports up to two dimensions)
%       .pl        path-loss factor
%       .c         speed of light in m/s
%       .f0        "center" frequency/frequencies for attenuation in Hz
%       .fs        sampling frequency in Hz
%
%   delay_s     delay in samples
%   gain        gain factor
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


function [delay_s, gain] = channel_large(settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   delay_s = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks / prepare parameters

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'type', 'dist', 'pl', 'c', 'f0', 'fs'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% settings.f0 might be a cell array => convert to matrix
if iscell(settings.f0)
   settings.f0 = cell2mat(settings.f0);
end

% check size of settings.dist and settings.f0
% (settings.f0 might contain entries of "deactivated" inputs and can thus be larger than settings.dist)
if length(settings.f0) < size(settings.dist,2)
   err('Not enough entries in settings.f0 or settings.dist too small.');
end

% match settings.f0 [in] to settings.dist [in, out]
settings.f0 = repmat(settings.f0(1:size(settings.dist,2)), size(settings.dist,1), 1);


% *******************************************************************************************************
% calculate attenuation and delay

% attenuation
switch lower(settings.type)
   case 'log-dist' % log-distance model
      gain =  ( settings.c ./ (4*pi* settings.dist .* settings.f0) ) .^ (settings.pl/2);
%    case 'log-shad' % log-distance model plus log-normal shading
   otherwise
      err('Unsupported large-scale type setting settings.type=''%s''', settings.type);
end
%     saturate to 1 (no active channels)
if any(any(gain > 1))
   critwarn('Found channel gains > 1 (active channel). Saturating gains to one.');
   gain(gain > 1) = 1;
end

% delay
delay_s = round(settings.dist / settings.c * settings.fs ); % samples
