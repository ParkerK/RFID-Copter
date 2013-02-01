% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% tag - demodulation plus sampling
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
% output = tag_demodulation(signal, clock, settings)
%    Models the demodulator as defined in demodulator.pdf (2007-08-27) and demodulator_vrfpsweep.doc
%    (2007-10-15). Demodulates signal according to settings and returns the sampled and demodulated
%    (binary) signal in output.
%
%
% ***** Interface definition *****
% function output = tag_demodulation(signal, clock, settings)
%    signal     decoder input signal (i.e. the modulated received signal)
%    clock      tag clock intervals (at sampling frequency) as created by TAG_CLOCK in mode 'fs'
%    settings   struct containing the demodulator settings
%       .fs               sampling frequency in Hz
%       .fc               carrier (center) frequency in Hz
%       .fclk             tag clock (center) frequency in Hz
%       .agc_len          periods of fc window length for AGC power detection (fcn. EST_POWER)
%       .agc_ol           overlapping in percent between AGC power detection windows (fcn. EST_POWER)
%       .agc_nbins        number of bins for AGC amplitude histogram
%       .tau              time constant of input RC filter in s
%       .force_exact_rc   force recalculation of RC filter coeff. if loaded ones are not an exact match
%       .dvpeak           steady state venv - vpeak
%       .ar               rise factor for vpeak (empirical)
%       .af               fall factor for vpeak
%       .slewrate         maximum slew rate for vpeak in V/s
%       .vdda             power supply (Vdda) voltage (scalar); Vgnd = Agnd = 0 V
%       .agc_sat          if true, out of range values of Vdda and carrier power are saturated to last
%                         available levels of AGC characteristic => demodulator is kept operable
%       .debouncelength   length of debounce filter (MA filter)
%       .turbo            downsampling prior to vpeak calculation if true (use with care!)
%       .turbo_osfactor   oversampling rate (critical sampling: 2*fclk) for "turbo" downsampling
%
%   output   decoded signal (binary @ fclk)
%
%
% ***** Used characteristics files (lookup tables) *****
% tagchar_demodulator.mat   AGC steady-state output voltage venv vs. vrf (peak) and vdda in [V]
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
% ? ACG power detection for 2 equally strong carriers
% ? time-variant AGC (->power supply fluctuations, carrier fluct)
% ? interpolation for RC filter coefficients if sampling frequency is not an exact match
%
% *******************************************************************************************************

function output = tag_demodulation(signal, clock, settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   output = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% downsampling for speed reasons ("turbo"): thresholds for warnings
internalsettings.osr_min   =  8; % minimum oversampling rate
internalsettings.rcatt_min = 30; % dB minimum attenuation of RC lowpass at new fs/2


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% number of input parameters
if nargin < 2
   criterr('Not enough input arguments.');
end

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'fs', 'fc', 'fclk', 'agc_len', 'agc_ol', 'agc_nbins', 'tau', 'force_exact_rc',...
   'dvpeak', 'ar', 'af', 'slewrate', 'agc_sat', 'debouncelength', 'vdda', 'turbo', 'turbo_osfactor'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% empty signal
if isempty(signal) || var(signal) == 0
   critwarn('Length or variance of input signal is zero. Returning empty vector.');
   output = [];
   return;
end

% -> column vectors
signal = signal(:);


% *******************************************************************************************************
% load characteristics (lookup tables)

% load AGC amplitude characteristic
tagchar_demod = loadmat('tagchar_demodulator');


% *******************************************************************************************************
% demodulator system: env (rectifier, AGC, RC-lowpass)
%    plus downsampling for performance reasons (if activated)

% downsampling factor (for performance reasons)
settings.dsf = round(settings.fs / settings.fclk / (settings.turbo_osfactor*2));

% get carrier peak amplitude for AGC
%     setup peakampl
settings.peakrms.fs     =  settings.fs; % sampling frequency in MHz
settings.peakrms.dsmode =           ''; % downsampling mode: no downsampling
settings.peakrms.dsf    = settings.dsf; % downsampling factor
settings.peakrms.nwin   = round(settings.fs / settings.fc * settings.agc_len); % window length ...
settings.peakrms.ol     = ceil(settings.peakrms.nwin * settings.agc_ol/100) / settings.peakrms.nwin * 100; % ... and overlapping
%     estimate (sine: peak = sqrt(2) * rms))
carrier_peak = sqrt(2) * peakrms(signal, settings.peakrms);

% rectifier
signal(find(signal < 0)) = 0; %#ok<FNDSB>

% RC lowpass, AGC => env
%     find closest index in RC filter coefficients
tagchar_demod.pos_fs = interp1(tagchar_demod.rc_fs, 1:tagchar_demod.settings.n_rc, settings.fs);
tagchar_demod.ind_fs = round(tagchar_demod.pos_fs);
%     create new RC filter coefficients if necessary (settings.tau=RC), use provided if possible (and
%     recalculation is case of a slight mismatch is not forced)
if isnan(tagchar_demod.pos_fs) || settings.tau ~= tagchar_demod.settings.rc_tau ||...
      ( settings.force_exact_rc && (tagchar_demod.ind_fs ~= tagchar_demod.pos_fs) )
   % make sure that the necessary toolbox is available (used only here => may not be available)
   if ~exist('tf', 'builtin')
      if isnan(tagchar_demod.pos_fs) || settings.tau ~= tagchar_demod.settings.tau
         err('Provided RC filter coefficients not applicable but cannot recalculate (toolbox missing).')
      else
         err('Cannot recalculate RC filter coefficients (toolbox missing). Consider switching off forced recalc.')
      end
   end
   % use c2d to discretize filter and finally tfdata to get filter coefficients
   [rc_b, rc_a] = tfdata(c2d(tf(1, [settings.tau, 1]), 1/settings.fs), 'v');
else
   % warn if there is no exact match
   if ~isnan(tagchar_demod.pos_fs) && (tagchar_demod.ind_fs ~= tagchar_demod.pos_fs)
      warn('No exact match for RC filter coefficients. Taking closest match (df=%.1f MHz).',...
         abs(settings.fs-tagchar_demod.rc_fs(tagchar_demod.ind_fs)));
   end
   rc_a = tagchar_demod.rc_a(tagchar_demod.ind_fs, :);
   rc_b = tagchar_demod.rc_b(tagchar_demod.ind_fs, :);
end
%     filter
signal = filter(rc_b, rc_a, signal);

% downsampling for performance reasons
if settings.turbo
   % calculate attenuation of RC lowpass at new fs/2
   rc_nyquistatt    = -20*log10(abs(freqz(rc_b, rc_a, [0, settings.fs/(2*settings.dsf)], settings.fs)));
   rc_nyquistatt(1) = [];
   % issue a message/warning
   msg('Downsampling by factor %i (new fs = %g x fclk), att of RC-lowp at new fs/2: %.1f dB.',...
      settings.dsf, settings.turbo_osfactor*2, rc_nyquistatt);
   if rc_nyquistatt < internalsettings.rcatt_min || settings.turbo_osfactor < internalsettings.osr_min
      critwarn('Low attenuation or critical oversampling factor.');
   end
   % downsample
   signal = signal(1 : settings.dsf : end);
   % adjust sampling clock
   clock = round(clock / settings.dsf);
   % make sure clock can still be used as index (...rounding)
   if clock(1) == 0
      clock(1) = 1;
   end
   if clock(end) > length(signal)
      clock(end) = length(signal);
   end
   % adjust settings
   settings.fs = settings.fs / settings.dsf;
   settings.ar = settings.ar * settings.dsf;
   settings.af = settings.af * settings.dsf;
end

% AGC (ideal)
%     warnings for out-of-range values plus saturation for "forced" operation
if settings.vdda < tagchar_demod.vdda(1) || settings.vdda > tagchar_demod.vdda(end)
   % forced operability: keep demodulator functional under all circumstances
   if settings.agc_sat
      critwarn('Vdda (%.2gV) outside characteristic (%.g-%.gV). Saturating.',...
         settings.vdda, tagchar_demod.vdda(1), tagchar_demod.vdda(end));
      settings.vdda = sanitize_setting('', settings.vdda, [tagchar_demod.vdda(1), tagchar_demod.vdda(end)]);
   % standard operations
   else
      warn('Vdda (%.2gV) outside characteristic (%.g-%.gV). Returning empty vector.',...
         settings.vdda, tagchar_demod.vdda(1), tagchar_demod.vdda(end));
      output = []; return;
   end
end
if carrier_peak < tagchar_demod.vrf(1) || carrier_peak > tagchar_demod.vrf(end)
   % forced operability: keep demodulator functional under all circumstances
   if settings.agc_sat
      warn('Signal peak amplitude (%.4gV) outside characteristic (%.g-%.gV). Saturating.',...
         carrier_peak, tagchar_demod.vrf(1), tagchar_demod.vrf(end));
      carrier_peak = sanitize_setting('', carrier_peak, [tagchar_demod.vrf(1), tagchar_demod.vrf(end)]);
      % standard operations
   else
      warn('Signal peak amplitude (%.4gV) outside characteristic (%.g-%.gV). Returning empty vector.',...
         carrier_peak, tagchar_demod.vrf(1), tagchar_demod.vrf(end));
      output = []; return;
   end
end
%     interpolate
[x, y] = meshgrid(tagchar_demod.vdda, tagchar_demod.vrf);
settings.agc_peak = interp2(x, y, tagchar_demod.venv, settings.vdda, carrier_peak);

% get signal peak amplitude (also for AGC)
%     setup peakampl
settings.peakrms.fs     =  settings.fs; % sampling frequency in MHz
settings.peakrms.dsmode =       'none'; % downsampling mode: random downsampling
settings.peakrms.nwin   = round(settings.fs / settings.fc * settings.agc_len); % window length ...
settings.peakrms.ol     = ceil(settings.peakrms.nwin * settings.agc_ol/100) / settings.peakrms.nwin * 100; % ... and overlapping
%     estimate (demodulated => ~DC)
demodulated_peak = peakrms(signal, settings.peakrms);

% scale to peak amplitude settings.agc_peak
signal = signal / demodulated_peak * settings.agc_peak;


% *******************************************************************************************************
% demodulator system: vpeak, output sampling by fclk and debouncing

% prepare parameters
settings.slewrate = settings.slewrate / settings.fs; % V/s / Hz => V/sample
level = settings.agc_peak - settings.dvpeak; % (assume carrier for t<0 => steady state already reached)

% generation of vpeak
vpeak = zeros(size(signal));
for i = 1 : length(signal) 
   % step to new level
   deltalevel = signal(i)-level-settings.dvpeak;
   % rise / fall speed
   if deltalevel > 0
      deltalevel = deltalevel * settings.ar;
   else
      deltalevel = deltalevel * settings.af;
   end
   % bounded slew rate
   if abs(deltalevel) > abs(settings.slewrate)
         deltalevel = sign(deltalevel) * settings.slewrate;
   end
   % new level
   level = level + deltalevel;
   % saturate
   if level > settings.vdda
      level = settings.vdda;
   end
   if level < 0 % vgnd
      level = 0;
   end
   % and store
   vpeak(i) = level;
end

% comparator and sampling
output = double( signal(clock) > vpeak(clock) );

% debouncing (filter shifts all edges, so this should not be a problem for timing detection)
debouncefilter = ones(settings.debouncelength,1) / settings.debouncelength;
output = round(filter(debouncefilter, 1, output));


return
% *******************************************************************************************************
% DEBUG (used to set time constants etc.)

% % do downsampling now (for plots) if not done before
% if ~settings.turbo
%    % downsample
%    signal = signal(1 : settings.dsf : end);
%    % adjust sampling clock
%    warn('Attributes of tag_clock changed by downsampling!'); % test makes not much sense => not done
%    clock = round(clock / settings.dsf);
%    % adjust settings
%    settings.fs = settings.fs / settings.dsf;
% end

% % signals plot
% t = (0:1:length(signal)-1) / settings.fs * 1e3;
% figure; hold on;
% plot(t, signal, 'k')
% plot(t, ones(size(signal))*settings.agc_peak, 'b--');
% plot(t, vpeak, 'r')
% plot(t(clock), output, 'g')
% hold off; grid minor;
% legend('env', 'desired env peak', 'vpeak', 'out');
% xlabel('t [ms]');

% % AGC (carrier amplitude) histogram
% figure;
% bar(ampl_bins, ampl_hist);
% title('AGC input amplitude histogram');
% xlabel('vrf [V]');

% % AGC (carrier amplitude) plot
% figure; hold on;
% t = (0:1:length(carrier_ampl)-1) / settings.fs;
% plot(t, carrier_ampl, 'b-');
% plot(t, ones(size(carrier_ampl)) * carrier_peak, 'k-');
% plot(t, ones(size(carrier_ampl)) * settings.agc_peak, 'r-');
% hold off; grid on;
% title('AGC functionality')
% legend('carrier peak amplitude', 'detected carrier peak level (keep limited histogram resolution in mind!)', 'AGC output level')
% xlabel('t [us]');

% AGC characteristic plot
% figure; hold on;
% surface(tagchar_demod.vrf, tagchar_demod.vdda, tagchar_demod.venv);
% hold off; grid on; view([-60, 15]);
% title('AGC characteristic (AGC output voltage)'); xlabel('vrf [V] (peak)'); xlabel('vdda [V]'); zlabel('venv [V] (peak)');
% colorbar();

% RC filter frequency response
%    continuous-time filter
% figure; bode(tf(1, [settings.tau, 1])); grid on;
%     discrete-time filter
% [rc_h, rc_f] = freqz(rc_b, rc_a, linspace(0, settings.fs/(2*settings.dsf), 1e3), settings.fs); 
% figure; plot(rc_f, 20*log10(abs(rc_h))); grid on; 
% set(gca, 'xScale','log'); xlim([min(rc_f), max(rc_f)]); ylim(xyzlimits(20*log10(abs(rc_h))))
% setlabels('RC LOWPASS FREQUENCY RESPONSE', 'f [Hz]   -   far right: new fs/2', 'H_{RC} [dB]');

