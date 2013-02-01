% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - time-variant (windowed) parameter estimation of a sinusoid (COSINE)
%         cf. Kay, Statistical Signal Processing - Estimation Theory, p. 168
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
% version = est_sinusoid()
%    Just returns the version number (string).
% param = est_sinusoid(signal, settings)
%    Estimates magnitude and phase of a presumed cosine (frequency SETTINGS.F0) in SIGNAL and returns
%    these parameters in PARAM.
%
%
% ***** Interface definition *****
% function param = est_sinusoid(signal, settings)
%    signal     signal vector
%    settings   struct containing settings
%       .n         length of windows in samples
%       .ol        overlapping of windows in percent
%       .f0        normalized frequency of sinusoid (fsin/fs)
%
%    param      struct containing the estimated parameters
%       .a         amplitude
%       .b         phase in degree (ATTENTION: PHASE RELATIVE TO A COSINE)
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
% ? add variance of estimates 
% ? other window that rectangular
%
% *******************************************************************************************************

function param = est_sinusoid(signal, settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   param = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'n', 'ol', 'f0'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end


% *******************************************************************************************************
% partitioning and calculation

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
settings.n_nol = floor(settings.n * (1 - settings.ol/100)); % better more overlapping than less
%     check for 100% overlapping (=> infinite number of windows)
if settings.n_nol == 0
   warn('Settings lead to 100%% overlapping. Setting non-overlapping part to 1 sample.');
   settings.n_nol = 1;
end
settings.n_ol  = settings.n - settings.n_nol;

% number of windows
settings.nwin = floor( (length(signal) - settings.n) / settings.n_nol ) + 1;  

% basis functions
cosine = cos(2*pi*settings.f0*[0:1:settings.n-1]');
sine   = sin(2*pi*settings.f0*[0:1:settings.n-1]');

% estimators
%     prepare array [alpha1, alpha2]
alpha = zeros(settings.nwin, 2);
%     estimators
for i = 0 : settings.nwin - 1
   alpha(1+i, 1) = 2 / settings.n * sum( signal(i*settings.n_nol+1:i*settings.n_nol+settings.n) .* cosine);
   alpha(1+i, 2) = 2 / settings.n * sum( signal(i*settings.n_nol+1:i*settings.n_nol+settings.n) .* sine  );
end
%     phase / amplitude
param.a = sqrt( alpha(:,1).^2 + alpha(:,2).^2 );
% param.p = - 180/pi * atan( alpha(:,2) ./ alpha(:,1) );
param.p = - 180/pi * atan2(alpha(:,2), alpha(:,1));
