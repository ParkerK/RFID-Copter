% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% tag - clock (sample time vector generator)
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
% version = tag_clock()
%    Just returns the version number (string).
% nt = tag_clock(settings)
%    Returns a vector containing sampling indices (nT + jitter).
%    Thus sampling can for example be done by sampled_vector = original_vector(tag_clock(settings)) 
%    with settings.length = length(original_vector) and settings.mode = 'fs'.
%
%
% ***** Interface definition *****
% function nt = tag_clock(settings)
%    settings     struct containing clock parameters
%       .fcenter     center ("clock") frequency in Hz
%       .fsigma      sampling time variation (sigma) in percent
%       .phi0        phase shift in deg (0<=phi0<360) ; -1: random
%       .fs          sampling frequency in Hz
%       .length      maximum length in samples (@fs or @fclk depending on mode)
%       .mode        {'fs', 'fclk'}
%                       'fs':   nt(end) <= settings.length (for downsampling)
%                       'fclk': length(nt) == settings.length
%
%    nt           vector containing sample indices
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


function nt = tag_clock(settings)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   nt = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'fcenter', 'fsigma', 'phi0', 'fs', 'length', 'mode'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% settings.length natural number and > 0 ?
if round(settings.length) ~= settings.length || settings.length < 0
   criterr('settings.length has to be a natural number > 0');
end

% settings.length == 0
if settings.length == 0
   warn('A clock vector of length zero was requested. Returning an empty vector.')
   nt = [];
   return
end


% *******************************************************************************************************
% create sampling intervals (samples @ fs)

% normalized times
tcenter = settings.fs / settings.fcenter;
tsigma  = tcenter * settings.fsigma / 100; % (settings.fsigma is sigma in percent nominal)

% phase shift -> delay
if settings.phi0 < 0 % if random
   settings.phi0 = 360 * rand;
end
settings.phi0 = mod(settings.phi0, 360); % truncate to one period
t0 = settings.phi0 / 360 * tcenter;

% create array
switch lower(settings.mode)
   case 'fs'
      % first and last entry: no jitter to avoid nt < 0 and nt > settings.nmax (at least for sane tsigma)
      length = floor(settings.length/tcenter);
      nt = 1 + round( t0 + (0:1:length-1)'*tcenter + [0; randn(length-2,1); 0] * tsigma );
   case 'fclk'
      % first entry: no jitter to avoid nt < 0 (at least for sane tsigma)
      length = settings.length;
      nt = 1 + round( t0 + (0:1:length-1)'*tcenter + [0; randn(length-1,1)] * tsigma );
   otherwise
      err('Unknown settings.mode = "%s"', settings.mode);
end

% check first and last entry (last entry only for mode 'fs': max(nt) <= settings.length)
if min(nt) < 1 || ( strcmp(settings.mode, 'fs') && max(nt) > settings.length )
   err('Generated sampling intervals (indices) contain entries( settings.nmax < nt(i) < 1). Check settings.');
end

% monotony check
if min(diff(nt)) <= 0
   critwarn('Generated sampling intervals (indices) are not monotonous. Check settings.');
end
