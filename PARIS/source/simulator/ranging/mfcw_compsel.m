% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% ranging - MFCW component selection according to the complex MF-CW derivation
%
% Note: Assumes \omega_0 = 0, hence SETTINGS.FI(i) = \omega_i, but c_mi(1) corresponds to \omega_0!
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
% version = mfcw_compsel()
%    Just returns the version number (string).
% Hd = mfcw_compsel(settings)
%     Returns the filter structure of the selection filter (lowpass).
%     IMPORTANT: filtfilt used => |Hd| has to be squared, arg(Hd) can be ignored (zero-phase)
% [c_mi, c_i, c_im] = mfcw_compsel(signal, settings)
%    Selects components according to the MF-CW derivation using IIR filters and returns their complex
%    amplitudes. Mixing is done with zero phase shift, all filters have zero group delay.
%   
%
%
% ***** Interface definition *****
% function [c_mi, c_i, c_im] = mfcw_compsel(varargin)
%    sig_bb     demodulated complex baseband signal
%    settings   struct containing settings
%       .nc        number of secondary carriers (has to be <= length(.fi))
%       .fi        array with secondary carrier frequency offsets (to main carrier: f0) in Hz
%       .lf        backscatter link frequency in Hz
%       .nharm     number of harmonics to consider for overlapping harmonics check 
%                  (non-sinusoidal tag modulation)
%       .iirord    filter order of selection filters
%       .att       stopband attenuation of selection filter in dB
%       .bw        stopband edge frequency @.fi+/-.bw in Hz
%       .frs       (reader) sampling frequency in Hz
%
%    c_mi    signal component amplitudes at omega_i+omega_m (assuming omega_0=0=>DC)
%            -> alternatively: filter object (overall dfilt filter structure) 
%    c_i     signal component amplitudes at omega_i (assuming omega_0=0=>DC)
%    c_im    signal component amplitudes at omega_i+omega_m (assuming omega_0=0=>DC)
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
% - santitize settings to prevent overlapping harmonics in the first place
%
% *******************************************************************************************************


function [c_mi, c_i, c_im] = mfcw_compsel(varargin)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   c_mi = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% demux input parameters
if nargin == 1 % return filter object
   signal   = [];
   settings = varargin{1};
   settings.justgauge = true;
elseif nargin == 2 % normal operation
   signal   = varargin{1}(:);
   settings = varargin{2};
   settings.justgauge = false;
else
   err('Wrong number of input arguments.');
end
   
% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'nc', 'fi', 'lf', 'nharm', 'iirord', 'att', 'bw', 'frs'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% no cw-ranging information whatsoever...
if settings.nc == 0
   msg('No MFCW secondary carriers. Is this a mistake?');
end

% check number of secondary carriers to add
if settings.nc > length(settings.fi) || settings.nc < 0
   err('Cannot add settings.nc=%i carriers; only %i defined in settings.fi and settings.vari.',...
      settings.nc, length(settings.fi));
end

% check frequency setup
if max(settings.fi) + settings.lf >= settings.frs/2
   err('Maximum cutoff (fi + fm) larger than frs/2. Check frequency setup.')
end

% empty signal
if isempty(signal)
   warn('Signal is empty. Returning NaNs.');
   c_mi = nan(1, settings.nc+1);
   c_i  = nan(1, settings.nc+1);
   c_im = nan(1, settings.nc+1); 
   return;
end
   
% use signal as leadin and leadout (filtfilt initialization)
signal = [signal(end:-1:1); signal; signal(end:-1:1)];


% *******************************************************************************************************
% check if harmonics overlap ... this will find all overlappings (also between nonadjacent carriers)

% calculate all harmonics (including carrier f_i=0)
settings.fi = [0; settings.fi(:)];
harmonics = [];
for i = 1 : settings.nc+1
      % harmonics for this frequency
      hi = abs(settings.fi(i))+settings.lf*[-settings.nharm:1:-1,1:1:settings.nharm]';
      % add only freq > 0 (symmetrical spectrum)
      harmonics = [harmonics; hi(hi>0)];
end
settings.fi(1) = [];

% sort and check smallest distance
harmonics = sort(harmonics);
if min(diff(harmonics)) < settings.bw
   critwarn('Overlapping harmonics possible (diff = %.1f kHz); filter might not be effective.', min(diff(harmonics))*1e-3);
end


% *******************************************************************************************************
% select components according to complex MF-CW derivation
% ... do not treat w_0 seperately (for simplicity, although likely settings.nc<<10)

% create DC filter (will be reused)
[b,a] = cheby2(settings.iirord, settings.att/2, settings.bw*2/settings.frs);
% figure; freqz(b,a, linspace(0,settings.lf,1e3), settings.frs); title('DC BANDPASS (ATTENUATION x2 BECAUSE OF filtfilt!)');
%     return filter object?
if settings.justgauge
   c_mi = dfilt.df1t(b, a);
   return;
end

% create mixer signals
%     time vector for mixers, correct phase shift caused by leadin (t=0 @ original sig_i,sig_q)
t = [-length(signal)/3:2/3*length(signal)-1]' / settings.frs; % len(sig_i) is 3 times length of before => /3 no problem
%     add w_0 = 0 to settings.fi (MFCW_ADDSECCARRIERS assumes w_0=0 and thus does not add this carrier)
settings.fi = [0; settings.fi(:)];

% shift to DC and filter (i=1 corresponds to \omega_0=0 here => around DC) 
for i = 1 : settings.nc + 1 % ... hopefully faster than vectorization, certainly less RAM
   c_mi(:,i) = filtfilt( b,a, signal .* exp(complex(0,2*pi*(-settings.fi(i)+settings.lf)*t)) ); % w_i - w_m
   c_i (:,i) = filtfilt( b,a, signal .* exp(complex(0,2*pi*(-settings.fi(i)            )*t)) ); % w_i
   c_im(:,i) = filtfilt( b,a, signal .* exp(complex(0,2*pi*(-settings.fi(i)-settings.lf)*t)) ); % w_i + w_m
end

% remove leadin/leadout
c_mi = c_mi(end/3+1:2*end/3, :);
c_i  = c_i (end/3+1:2*end/3, :);
c_im = c_im(end/3+1:2*end/3, :);
