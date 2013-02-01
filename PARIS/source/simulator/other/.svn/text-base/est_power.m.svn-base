% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - time-variant (windowed) power (mean^2 + variance) calculation of signals
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
% version = est_power()
%    Just returns the version number (string).
% power = est_power(signal, settings)
%    Calculated the time-variant power of SIGNAL according to SETTINGS (variance of SIGNAL within 
%    rectangular windows of length SETTINGS.N, overlapping by SETTINGS.OL percent).
%
%
% ***** Interface definition *****
% function power = est_power(signal, settings)
%    signal     signal vector
%    settings   struct containing settings
%       .n         length of windows in samples
%       .ol        overlapping of windows in percent
%
%    power   vector containing the (time-variant) variance
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
% - merge with est_sinusoid to a generic windowed estimator (nearly idential function)
%   + add variance of estimate
% ? other window that rectangular
%
% *******************************************************************************************************

function power = est_power(signal, settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   power = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'n', 'ol'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end


% *******************************************************************************************************
% partitioning and variance calculation

% check if settings.n is integer
if settings.n ~= round(settings.n)
   warn('Window length (settings.n) is not integer: rounding.');
   settings.n = round(settings.n);
end

% window longer than signal?
if settings.n > length(signal)
   warn('Window length (settings.n) is larger than provided signal. Setting window size to signal length.');
   settings.n = length(signal);
end

% length of (non-)overlapping part of window
settings.n_nol = round(settings.n * (1 - settings.ol/100)); % better more overlapping than less
%     check for 100% overlapping (=> infinite number of windows)
if settings.n_nol == 0
   warn('Settings lead to 100%% overlapping. Setting non-overlapping part to 1 sample.');
   settings.n_nol = 1;
end
settings.n_ol  = settings.n - settings.n_nol;

% number of windows
settings.nwin = floor( (length(signal) - settings.n) / settings.n_nol ) + 1;  

% calculate variance
power = zeros(settings.nwin, 1);
for i = 1 : settings.nwin
   power(i) = mean(signal( (i-1)*settings.n_nol+1 : (i-1)*settings.n_nol+settings.n))^2 + ... % DC
      var( signal( (i-1)*settings.n_nol+1 : (i-1)*settings.n_nol+settings.n) ); % AC
end

