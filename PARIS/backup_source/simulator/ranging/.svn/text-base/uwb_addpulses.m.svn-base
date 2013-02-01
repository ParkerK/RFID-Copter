% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% ranging - UWB pulse sequence generation
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
% version = uwb_addpulses()
%    Just returns the version number (string).
% [signal, pulse] = uwb_addpulses(signal, settings, oscsettings)
%
%
% ***** Interface definition *****
% function signal = mfcw_addseccarriers(signal, settings)
%    signal     main carrier signal
%    settings   struct containing settings
%       .fs        sampling frequency in Hz
%       .f0        center frequency in Hz
%       .bw        bandwidth of each pulse in Hz
%       .tpr       pulse repetition time in s
%       .shape     shape of added pulses {'rrc': root-raised-cosine}
%       .r         parameter for the pulses (rolloff factor 0..1 for shape 'rrc')
%       .var       pulse variance (relative to variance of SIGNAL)
%       .range_s   areas of signal to add secondary carriers to; leave empty to select entire SIGNAL
%                  (in samples @ fs [start1, end1, start2, end2, ...])
%    oscsettings   struct containing settings for READER_OSCILLATOR (modulation signal)
%
%   signal   main carrier signal with added UWB pulses
%   pulse    a single pulse (real baseband; can be used for matched filtering)
%   pos      pulse positions in samples @ fs
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
% - ternary pulse sequence (preamble)
% ? bandwidth correct
%
% *******************************************************************************************************


function [signal, pulse, pos] = uwb_addpulses(signal, settings, oscsettings)
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
% internal settings

% length of pulse templates
internalsettings.nper = 10; % times sampling frequency divided by bandwidth


% *******************************************************************************************************
% input parameter checks

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'fs', 'f0', 'bw', 'tpr', 'shape', 'r', 'var', 'range_s'};
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
% add a UWB pulse sequence to signal

% generate a single pulse
%     shape
switch lower(settings.shape)
   case 'rrc' % root raised cosine
      settings.df = settings.bw/2 * settings.r;
      pulse = firrcos(round(internalsettings.nper*settings.fs/settings.bw), settings.bw/2-settings.df/2, settings.df, settings.fs, 'sqrt');
   otherwise
      err('Undefined pulse shape "%s".', lower(settings.shape));
end
%     normalize variance
pulse = pulse / std(pulse) * settings.var * std(signal);

% generate pulse sequence
%     ideal
pulse_seq = zeros(size(signal));
pos = round(1 : settings.fs*settings.tpr : length(pulse_seq)); % pulse positions
pulse_seq(pos) = 1; % place pulses
%     apply pulse shape
pulse_seq = filter(pulse, 1, pulse_seq);

% modulate (shift in frequency) and add to signal
%    ... this is done for the pulse sequence instead of the pulse to ensure a continuous phase 
%        for the modulation carrier
%     create indices for UWB pulse ranges (length of range_s is even)
ind = [];
for i = 1 : length(settings.range_s)/2
   ind = [ind; [settings.range_s(2*i-1):settings.range_s(2*i)]'];
end
%     modulate/add
oscsettings.length  = length(signal) / settings.fs;
modcarrier = reader_oscillator(oscsettings);
signal(ind) = signal(ind) + pulse_seq(ind) .* modcarrier(ind);

% remove pulse positions that are outside settings.range_s
for i =  1 : length(pos)
   if ~any(ind == pos(i))
      pos(i) = NaN;
   end
end
pos(isnan(pos)) = [];
