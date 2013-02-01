% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% ranging - MFCW secondary carrier generation according to complex MFCW derivation
%
% Note: Assumes \omega_0 = 0, hence SETTINGS.FI(i) = \omega_i
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
% version = mfcw_addseccarriers()
%    Just returns the version number (string).
% signal = = mfcw_addseccarriers(signal, settings)
%    Secondary carriers according to SETTINGS are added to SIGNAL. Note that the variance of SIGNAL has
%    to be nonzero. The phase of each secondary carrier is zero for signal(1) and furthermore does not
%    restart at zero for each new block in SETTINGS.RANGE_S.
%   
%
%
% ***** Interface definition *****
% function signal = mfcw_addseccarriers(signal, settings)
%    signal     main carrier signal
%    settings   struct containing settings
%       .fs        sampling frequency in Hz
%       .f0        carrier frequency in Hz
%       .fi        array with secondary carrier frequency offsets (to main carrier: f0) in Hz
%       .vari      array with secondary carrier variances (relative to variance of SIGNAL)
%       .gammai    phase shift of secondary carriers [rad]
%       .nc        number of secondary carriers to add (has to be <= length(.fi)==length(.vari))
%       .range_s   areas of signal to add secondary carriers to; leave empty to select entire SIGNAL
%                  (in samples @ fs [start1, end1, start2, end2, ...])
%
%   signal   main carrier signal with added secondary carriers
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


function signal = mfcw_addseccarriers(signal, settings)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   signal = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'fs', 'f0', 'fi', 'vari', 'nc', 'range_s'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% add nothing?
if settings.nc==0 || isempty(settings.fi) || isempty(settings.vari)
   warn('mfcw_addseccarriers has been called with zero carriers to add. Is this a mistake?');
   return
end

% no main carrier?
if var(signal) == 0
   err('Main carrier (input parameter "signal") has zero variance.');
end

% range
%     empty => entire signal
if isempty(settings.range_s)
   settings.range_s = [1,length(signal)];
end
%     has to be inside signal
if any(settings.range_s) < 1 || any(settings.range_s > length(signal)) 
   err('Given settings.range_s exceeds length of signal.');
end
%     monotony
if any(diff(settings.range_s) < 1)
   err('settings.range_s is not monotonous.');
end
%     length has to be even [start1, end1, start2, end2, ...]
if mod(length(settings.range_s), 2) ~= 0
   err('Length of settings.range_s has to be even [start1, end1, start2, end2, ...].');
end

% other settings, lengths
if length(settings.fi) ~= length(settings.vari)
   err('Size of settings.fi and settings.vari does not match.');
end
if settings.nc > length(settings.fi) || settings.nc < 0
   err('Cannot add settings.nc=%i carriers; only %i defined in settings.fi and settings.vari.',...
      settings.nc, length(settings.fi));
end


% *******************************************************************************************************
% add secondary carriers to signal

% make sure signal is a column vector and calculate variance
signal = signal(:);
settings.var0 = var(signal);

% create indices (length of range_s is even)
ind = [];
for i = 1 : length(settings.range_s)/2
   ind = [ind; [settings.range_s(2*i-1):settings.range_s(2*i)]'];
end
t = (ind-1) / settings.fs; % [s], t=0 @ signal(1)

% add secondary carriers assuming perfect modulation (linear, noiseless)
% ... secondary carriers will likely be created by the digital part and be modulated onto the carrier
%     => perfect amplitude/frequency/phase (zero phase shift)
for i = 1 : settings.nc   
   signal(ind) = signal(ind) + cos(2*pi*(settings.f0+settings.fi(i)) * t) * sqrt(2*settings.vari(i)/settings.var0);
end
   
