% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% ranging - FMCW secondary carrier generation
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
% version = fmcw_addfmcarrier()
%    Just returns the version number (string).
% signal = fmcw_addfmcarrier(signal, settings)
%    Adds a frequency modulated secondary carrier to SIGNAL, according to SETTINGS. Note that the
%    variance of SIGNAL has to be nonzero. 
%   
%
%
% ***** Interface definition *****
% function signal = fmcw_addfmcarrier(signal, settings)
%    signal     main carrier signal
%    settings   struct containing settings
%       .fcn       frequency modulation function {'square', 'sawtooth', 'triangle', 'sin'}
%       .fs        sampling frequency in Hz
%       .f0        FM carrier center frequency in Hz
%       .df        maximum deviation from center frequency (+/-) during FM in Hz 
%       .per       period of one modulation cycle in s
%       .var       secondary carrier variance relative to signal variance
%       .range_s   areas of signal to add secondary carriers to; leave empty to select entire SIGNAL
%                  (in samples @ fs [start1, end1, start2, end2, ...])
%       .rfmod     RF carrier generation / modulation method
%                     'ideal':  generate the RF carrier directly (do not use the carrier in SIGNAL)
%                     'hilb':   single-sideband using a Hilbert transform to suppress the lower sideband
%
%   signal   main carrier signal with added secondary (FM-)carrier
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
% - multi-frequency-modulation CW ranging?
% - add series expansion or maximum significant modulation frequency / FM component
%
% *******************************************************************************************************


function signal = fmcw_addfmcarrier(signal, settings)
version = 'UC';

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
expected.reqfields = {'fcn', 'fs', 'f0', 'df', 'per', 'var', 'range_s', 'rfmod'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
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


% *******************************************************************************************************
% add frequency modulated secondary carrier to signal

% make sure signal is a column vector and calculate its variance
signal = signal(:);
settings.var0 = var(signal);

% create indices and time vector (length of range_s is even)
ind = [];
for i = 1 : length(settings.range_s)/2
   ind = [ind; [settings.range_s(2*i-1):settings.range_s(2*i)]'];
end
t = (ind-1) / settings.fs; % [s], t=0 @ signal(1)

% select modulation signal and calculate fourier series components [!!!] TO BE IMPLEMENTED
switch lower(settings.fcn)
   case 'square'
      xm =   square(2*pi/settings.per * t + pi);
   case 'sawtooth'
      xm = sawtooth(2*pi/settings.per * t + pi);
   case 'triangle'
      xm = sawtooth(2*pi/settings.per * t, 0.5);
   case 'sin'
      xm =      sin(2*pi/settings.per * t); 
   otherwise
      err('Unsupported modulation shape: "%s"', lower(settings.fcn));
end

% add FM secondary carrier using a rectangular approximation for the integral of xm
% ... will likely be created by the digital part and be modulated onto the carrier
%     => perfect and noiseless FM modulation (amplitude/frequency/phase)
switch(lower(settings.rfmod))
   
   % ideal/perfect RF carrier modulation: generate directly at settings.f0
   case 'ideal'
      signal(ind) = signal(ind) + sqrt(2*settings.var/settings.var0) * ...
         cos(2*pi*settings.f0*t + 2*pi*settings.df*cumsum(xm/settings.fs));

   % nonperfect RF carrier modulation: use carrier in signal for modulation and suppress LSB using hilb
   case 'hilb'
      fm = cos(2*pi*settings.df*cumsum(xm/settings.fs)) * sqrt(2*settings.var/settings.var0); % baseband
      signal(ind) = signal(ind) + signal(ind) .* fm + imag(hilbert(signal(ind))) .* imag(hilbert(fm));
      
   otherwise
      err('Unsupported method "%s".', settings.rfmod);
end
