% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% tag - modulation (inlcudes functionality of tag_power)
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
% version = tag_modulation()
%    Just returns the version number (string).
% settings = tag_modulation(data, settings)
%    Returns completed (e.g. the needed carrier length to encode and modulate DATA according to
%    SETTINGS).
% backscatter = tag_modulation(carrier, data, clock, settings)
%    Calculates the modulated backscatter signal according to SETTINGS. Transients are placed at the
%    clock ticks in CLOCK. No rolloff or additional filtering is done. If DATA contains only zeros
%    (CLOCK can be empty in this case), BACKSCATTER contains the unmodulated backscatter signal
%    (see next item).
% backscatter = tag_modulation(carrier, [0], [], settings)
%    Calculates the unmodulated backscatter signal according to SETTINGS. 
% transmitted = tag_modulation(carrier, [], [], settings)
%    Calculates the transmitted carrier according to SETTINGS. The returned vector is multiplied by the 
%    tag input impedance at the operating point (tag power, SETTINGS.FC => frequency flat!) and can  
%    be used as input to tag_demodulation directly ([Vrms] voltage drop over abs input impedance).
% backscatter = tag_modulation(carrier, modsignal, settings)
%    Calculates the modulated backscatter signal according to SETTINGS using an external modulation
%    signal. Used fclk=fs. Note that length calculations (e.g. needed carrier length) will not work
%    for this mode.
%
% [backscatter or transmitted, internal] = ...
%    Also returns interal values, for example, the estimated chip input power PIC.
%    This works for all modes above that return a carrier signal (not: version or settings).
%    
%
%
% ***** Interface definition *****
% function [modcarrier, internal] = tag_modulation(carrier, data, clock, settings)
%    carrier     (optional) input signal (= carrier) to modulate or filter
%    modsignal   (optional) external modulation signal (@fs); binary alphabet {0,1}
%    data        (optional) data to encode (bitstream); binary alphabet {0,1}
%    clock       (optional) tag clock intervals (at fs) as created by TAG_CLOCK in mode 'fclk'  
%    settings    struct containing the demodulator settings
%       .fs            sampling frequency in Hz
%       .fc            carrier (center) frequency in Hz
%       .fclk          tag clock (center) frequency in Hz
%       .lf            backscatter link frequency in Hz {40...640 kHz}
%       .t0            s of unmodulated carrier before and after modulation (not 100% exact!); 
%                      Note: .) t0 might be sanitized to avoid truncation of modulated carrier
%                            .) trailing space will be larger or equal settings.t0 if CLOCK has jitter
%       .charfile      filename for reflection coefficient characteristic .mat-file
%       .nfft          filterbank size (number of channels)
%       .nfilt         number of filterbank channel filter (FIR) coefficients
%       .imagtol       maximum average imag(backscatter.^2)/real(backscatter.^2) after filterbank synthesis 
%                      before a warning is issued (safety check; synthesized signal has to be real)
%       .force         forced modulation/filtering: if true, functionality checks are ignored (if possible) 
%       .ptag_f        force tag power level (set to NaN for normal operation)
%       .pfnc_bias     bias towards "unmodulated" for Pin and Pic (0...1)
%                      0: assume continuous operation (take power loss due to modulation into account)
%                      1: assume short-term operation (power supply buffer is able to deal w. mod. pwr loss)
%       .pwrcharfile   filename for power characteristic .mat-file
%       .pwrchar       power supply characteristic {'Sofie_V2A', 'Sofie_V2B', 'Sofie_V2C', 'Sofie_AVG'}
%                  
%   modcarrier   calculated backscatter signal or returned (sanitized, completed) settings (only if used
%                 for modulation, not filtering)
%      -> or: additional (sanitized) settings
%      .mode            detected mode {'modulate', 'filter'}
%      .t0              sanitized settings.t0 in s
%      .t0_s            sanitized settings.t0 in samples
%      .t_fclk_s        subcarrier period in samples @ fclk     
%      .length_fclk_s   length of modulated part in samples @ fclk (only modulation, no leadin/leadout)
%      .length_s        length of modulated part plus leadin/leadout in samples  ("needed carrier length")
%      .length          length of modulated part plus leadin/leadout in s ("needed carrier length")
%      .grpdel_s        group delay of filterbank (= removed number of samples)
%      .moddel_s        delay for modulation signal in samples (estimate; exact delay jitters)
%   internal     struct containing internal values, like, e.g., recorded power levels
%      .pav          available power in W
%      .pin          input power in W
%      .pic          chip input power in W
%      .vdda         power supply voltage in V 
%                    ATTENTION: will return the minimum Vdda for out-of-rangechip input powers Pic 
%      .fres         frequency resolution of filterbank in Hz
%      .maxerr_mod   max error caused by filterbank not including interpolation [abs, angle] for mod
%      .maxerr_umd   max error caused by filterbank not including interpolation [abs, angle] for unmod
%      .rho          center value of linear reflection coefficient model .rho(f)
%                    ... for frequencies f = linspace(0,settings.fs/2,settings.nfft/2+1)
%      .drho         difference value of linear reflection coefficient model .drho(f)
%                    ... for frequencies f = linspace(0,settings.fs/2,settings.nfft/2+1)
%      .zic          estimated chip input impedance not including modulation in Ohm
%                    (only for mode "filter", NaN otherwise) 
%      .zat          used assembly impedance in Ohm (only for mode "filter", NaN otherwise) 
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
% ? carrier peak power detection instead of average (or something in between)
%   (avg power leads to lower power for modulated signal => Zin>> => high output voltage for 'filter')
% ? random t0
% - store filter coefficients (optimized offline) and not reflection coefficient
%   + time-variant power detection
%   + smooth transients for Rmod
% - ultrawideband power detection / input impedance calculation for multiple carriers 
%   (remove dependence on settings.fc)
% = extract middle part to tag_power, add power supply buffer model
%
% *******************************************************************************************************


function [modcarrier, internal] = tag_modulation(varargin)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   modcarrier = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% filename of filterbank channel coefficients
internalsettings.fbcharfile = 'tagchar_modulator_fb';

% maximum filterbank length as fraction of min length("high", "low")
internalsettings.nfft_maxfraction = 5; % (< 2 is not a good idea)

% tolerances before a warning is issued
internalsettings.lf_tol  =  100; % Hz "fclk multiple of LF"
internalsettings.imagtol = 1e-4; % im/re after synthesis

% factor for maximum error created by finite filterbank resolution
% (approximate, cf. output of plots_tag_modulation.m) 
internalsettings.fbdiff = 1/3; % times diff(spectrum)


% *******************************************************************************************************
% input parameter checks

% number of input parameters
if nargin == 2
   justreturnsettings = true;
   data     = varargin{1}(:);
   settings = varargin{2};
elseif nargin == 3
   critwarn('External modulation signal. Setting fclk to fs.');
   justreturnsettings = false;
   carrier        = varargin{1}(:);
   modsignal_fclk = varargin{2}(:);
   settings       = varargin{3};
   settings.fclk  = settings.fs; % fclk==fs (fastest possible reaction)
   clock          = 1:length(modsignal_fclk);
   data           = [0]; % dummy
elseif nargin == 4
   justreturnsettings = false;
   carrier  = varargin{1}(:);
   data     = varargin{2}(:);
   clock    = varargin{3}(:);
   settings = varargin{4};
else
   criterr('Wrong amount of input arguments.');
end

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'fs', 'fc', 'fclk', 'lf', 't0', 'charfile', 'nfft', 'nfilt', 'force', 'ptag_f',...
   'pfnc_bias', 'pwrcharfile', 'pwrchar'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% length of carrier signal checked below (when timing is fixed)

% check if data signal has alphabet {0,1}
if sum(data == 0) + sum(data == 1) ~= length(data)
   err('Data has non-binary alphabet. Only binary signals with alphabet {0,1} allowed.')
end

% for external modulation signal: check if modsignal_fclk has alphabet {0,1}
if nargin == 3 && ( sum(modsignal_fclk == 0) + sum(modsignal_fclk == 1) ~= length(modsignal_fclk) )
   err('External modulation signal must have binary alphabet ({0,1}).')
end

% check settings.nfft (power of 2)
if log2(settings.nfft) ~= round(log2(settings.nfft))
   err('Filterbank FFT size (settings.nfft) has to be a power of 2.');
end

% mode?
%     filter (transmit signal)
if nargin == 4 && isempty(data) && isempty(clock)
   settings.mode = 'filter';
%     just reflect
elseif nargin == 4 && sum(data)==0 && isempty(clock)
   settings.mode = 'modulate';
   clock = [1, length(carrier)];
%     modulate data
else
   settings.mode = 'modulate';
end

% initialize struct internal
%     power levels
internal.pav = NaN;
internal.pin = NaN;
internal.pic = NaN;
%     power supply voltage
internal.vdda = NaN;
%     functional
internal.func =  0;
%     frequency resolution of filterbank
internal.fres = NaN;
%     maximum error caused by filterbank not including interpolation [abs, angle]
internal.maxerr_mod = [NaN, NaN];
internal.maxerr_umd = [NaN, NaN];
%     linear model
internal.rho  = NaN;
internal.drho = NaN;
%     impedances
internal.zic = NaN;
internal.zat = NaN;


% *******************************************************************************************************
% lengths (encoding takes care of data rates; here we are at subcarrier level; Smod = 0/1)

% modulation
if nargin ~= 3 && strcmpi(settings.mode, 'modulate')
   % check if data is empty (...likely a mistake)
   if isempty(data)
      critwarn('Data vector is empty.');
   end
   
   % at clock frequency fclk (rounding: fclk has to be high enough to satisfy LF tolerances) [???]
   %     issue warning if fclk is not a multiple of LF
   if strcmpi(settings.mode, 'modulate') && mod(settings.fclk, settings.lf) > internalsettings.lf_tol
      critwarn('Tag clock frequency is not a multiple of LF (original: %.1f kHz, corrected: %.1f kHz)', ...
         settings.lf, settings.fclk / round(settings.fclk/settings.lf));
   end
   %     subcarrier period
   settings.t_fclk_s = round(settings.fclk / settings.lf);
   %     length of modulated (0) and unmodulated (1) halfbit (make modulated halfbit always shorter -> energy) 
   length0_fclk_s = floor(settings.t_fclk_s / 2);
   length1_fclk_s = ceil(settings.t_fclk_s / 2);
   %     length of modulated part
   settings.length_fclk_s = (length(data)-sum(data)) * length0_fclk_s + sum(data)*length1_fclk_s;  % Smod==0 + Smod==1

   % at sampling frequency fs (rounding error here: don't care)
   %     length of modulated (0) and unmodulated (1) halfbits
   length0_s = round(length0_fclk_s * settings.fs / settings.fclk);
   length1_s = round(length1_fclk_s * settings.fs / settings.fclk);   
   
   % modulation with external modulation signal: no checks
elseif nargin == 3 && strcmpi(settings.mode, 'modulate')
   length0_s = length(carrier);
   length1_s = length(carrier);
   
% no modulation, just filter ... like modulating a zero that lasts for the whole carrier length
else
   % "the whole carrier lenght"
   clock = [1, length(carrier)];
   % before (and after) data, there are zeros
   data  = []; 
   % set t0 to zero (no unfiltered parts!)
   settings.t0 = 0;
   % no modulation patterns
   length0_fclk_s = 0;
   length1_fclk_s = 0;
   settings.length_fclk_s  = 0;
end

% load filterbank characteristics
tagchar_mod_filterbank = loadmat(internalsettings.fbcharfile);

% check if filter is longer than a fraction of backscatter link period
if strcmpi(settings.mode, 'modulate') && (settings.nfft > min(length0_s, length1_s)/internalsettings.nfft_maxfraction)
   nfft = 2^floor(log2(min(length0_s, length1_s)/internalsettings.nfft_maxfraction));
   critwarn('Filterbank too large for LF. Reducing channels to %i (old: %i).', nfft, settings.nfft);
   settings.nfft = nfft;
end

% length of "leadin" and "leadout"
%     length of t0 in samples
settings.t0_s = round( settings.t0 * settings.fs );
%     make sure modulated part is not truncated in case of clock phase shift and clock jitter
%     ... one tclk for phase shift, one tclk for jitter (on safe side) plus nfft/2
if strcmpi(settings.mode, 'modulate') && (settings.t0_s < 2*round(settings.fs / settings.fclk) + settings.nfft/2)
   critwarn('Saturating t0_s to minimum to prevent truncation of mod part.')
   settings.t0_s = 2*round(settings.fs / settings.fclk) + settings.nfft/2;
end
%     in seconds (... sanitized setting)
settings.t0 = settings.t0_s / settings.fs;

% group delay of filterbank
%     delay for carrier
settings.grpdel_s = settings.nfft/2 * ( settings.nfilt - 1);
%     delay for modulation signal
%        the modulation signal has a smaller delay than the carrier (grpdel_s); nontheless grpdel_s is
%        removed from the modulated signal => this has to be compensated by adding the delay difference
%         to the clock signal (plus block processing!)
%           ... modulation reaches ind_clk == 1 at round(t0_s/(nfft/2)) if set delay is settings.t0_s+
%               round(settings.fs/settings.fclk)-settings.nfft/4 iand initial clock jitter is 
%               compensated: -clock(1)+1
%            ... shift delay in that case is below settings.t0_s (removal of grpdel_s)
%            ... error to removed grpdel_s measured to 50% value of first transient is approximately
%                (settings.nfilt/2-2)*settings.nfft/2
%           => add delay for ind_clk==1@i*nfft/2==t0_s, add removed grpdel_s and compensate
%           ... grpdel_s - (nfilt/2-2)*nfft/2 = (nfilt/2+1)*nfft/2
%     ... problem: t0_s = 0 is dangerous (not "unmodulated" at beginning)
settings.moddel_s = round(settings.fs/settings.fclk) - settings.nfft/4 +... % ind_clk==1@i*nfft/2==t0_s
                    (settings.nfilt/2+1)*settings.nfft/2; % difference of real delay to rem. grpdel_s
                 
% just return settings if requested
if justreturnsettings
   % length of whole modulated part (estimate; will change because of clock jitter)
   settings.length_s = (length(data)-sum(data)) * length0_s + sum(data)*length1_s +... % Smod==0 + Smod==1
      2*settings.t0_s +... % additional (manual) leadin/leadout
      settings.grpdel_s; % group delay of filterbank (deleted part)

   %     make settings.lenths_s a multiple of settings.nfft/2 (for PPN filterbank)
   settings.length_s = settings.length_s + settings.nfft/2 - mod(settings.length_s, settings.nfft/2); % samples
   settings.length   = settings.length_s / settings.fs; % s

   % return these settings
   modcarrier = settings;
   return
end
  
% recalculate length @ fs (... clock jitter of tag clock)
% (Attention: length_s is length of modulated/filtered part and not carrier length!)
settings.length_s = clock(end); % samples
settings.length   = settings.length_s / settings.fs; % s

% check length of carrier signal now that timing is fixed (group delay is already included in settings.length_s)
%     make sure carrier length is a multiple of settings.nfft
if mod(length(carrier), settings.nfft/2) ~= 0
   critwarn('Length of carrier (%i) is no multiple of settings.nfft/2 (%i). Truncating.', length(carrier), settings.nfft/2);
   carrier(end-mod(length(carrier), settings.nfft/2)+1:end) = [];
end
%     check carrier length
if strcmpi(settings.mode, 'modulate') && (length(carrier) < settings.length_s)
   err('Insufficient length of carrier signal for given data vector (%i < %i).', length(carrier), settings.length_s);
end


% *******************************************************************************************************
% create modulation signal @ fclk (in range 0...1)

% build modulation signal at tag-clock level (only if generated internally)
if nargin ~= 3
   mod0_fclk = zeros(length0_fclk_s, 1);
   mod1_fclk = ones(length1_fclk_s, 1);
   modsignal_fclk = zeros(settings.length_fclk_s, 1);
   index = 1;
   %     create signal bit per bit
   for i = 1 : length(data)
      if data(i) == 0
         modsignal_fclk(index:index+length0_fclk_s-1) = mod0_fclk;
         index = index + length0_fclk_s;
      else
         modsignal_fclk(index:index+length1_fclk_s-1) = mod1_fclk;
         index = index + length1_fclk_s;
      end
   end
end


% *******************************************************************************************************
% characteristics: reflection coefficients of tag input + determine power levels and support

% load reflection coefficient characteristic
tagchar_mod = loadmat(settings.charfile);

% determine incident (available) carrier power level
if ~isnan(settings.ptag_f)
   tagchar_mod.pav_hat = settings.ptag_f;
   critwarn('Forced tag power level: %.2e uW (%.2f dBm).', tagchar_mod.pav_hat*1e6, 10*log10(tagchar_mod.pav_hat/1e-3));
else
   %     setup estimator (best estimation for carrier frequency)
   %     ... window/overlapping size multiples of carrier period 1/f0, ~50% overlapping
   %     ... window size max 1/4 t_lf samples (for time-variant: rho can be selected 4 times per t_lf)
   %         (if applicable)
   setup_est.f0 = settings.fc / settings.fs;
   setup_est.n  = round(ceil( settings.fc/(4*settings.lf) ) / setup_est.f0 );
   setup_est.ol = round(ceil( setup_est.n/2 * setup_est.f0 ) / setup_est.f0) / setup_est.n * 100;
   %     estimate (for now: only time-invariant, i.e. short-time stationary channel)
   %     ... average power => assumes that the current signal is representative and tag is well buffered
   tagchar_mod.pav_hat = mean(est_power(carrier, setup_est));
end
%     output
msg('Carrier power level (variance): Pav = %.2e uW (%.2f dBm)', tagchar_mod.pav_hat*1e6, 10*log10(tagchar_mod.pav_hat*1e3));

% get operating point (combination power/frequency: position, rounded index and vicinity)
%     available power
tagchar_mod.pos_pav = interp1(tagchar_mod.pav, [1:1:tagchar_mod.settings.np], tagchar_mod.pav_hat, 'linear',  NaN);
tagchar_mod.ind_pav = round(tagchar_mod.pos_pav);
if floor(tagchar_mod.pos_pav) == tagchar_mod.settings.np % very unlikely, but possible
   tagchar_mod.vic_pav = floor(tagchar_mod.pos_pav) + [0,0];
else % standard case
   tagchar_mod.vic_pav = floor(tagchar_mod.pos_pav) + [0,1];
end 
%     frequency (non-extrapolated)
tagchar_mod.pos_fch = interp1(tagchar_mod.fch, [1:1:tagchar_mod.settings.nfch], settings.fc, 'linear', NaN);
tagchar_mod.ind_fch = round(tagchar_mod.pos_fch);
if floor(tagchar_mod.pos_fch) == tagchar_mod.settings.nfch % very unlikely, but possible
   tagchar_mod.vic_fch = floor(tagchar_mod.pos_fch) + [0,0];
else % standard case
   tagchar_mod.vic_fch = floor(tagchar_mod.pos_fch) + [0,1];
end

% if frequency is outside characteristic: setup problem
if isnan(tagchar_mod.pos_fch)
   err('Carrier frequency (settings.fc) out of characteristic range (%.2f to %.2f MHz). Check setup.',...
   tagchar_mod.fch(1)/1e6, tagchar_mod.fch(end)/1e6);
end

% do we have characteristic support (if power in range of pav vector and no NaNs in vicinity of op.pt.) ?
%     extrapolated (if not: modulation/reflection will not work at all)
tagchar_mod.support.ext = ~isnan(tagchar_mod.pos_pav) && ...
   ~any(isnan(sum(sum(tagchar_mod.rhoext_pav(1+tagchar_mod.settings.nf_out1+tagchar_mod.vic_fch, :, tagchar_mod.vic_pav)))));
%     non-extrapolated (if not: tag is definitely not functional)
%     ... the nonextrap. char. contains values far below pmin => this cannot be covered by the pmin check
%         below
tagchar_mod.support.nonext = ~isnan(tagchar_mod.pos_pav) && ...
   ~any(isnan(sum(sum(tagchar_mod.rho_pav(tagchar_mod.vic_fch, :, tagchar_mod.vic_pav)))));

% calculate input power Pin and chip input power Pic
% ... note that tagchar_mod.p_in might be NaN for high power levels, even if there is support (interp1)
if tagchar_mod.support.ext % we need at least extrapolated support to do this
   % we are doing a square modulation (below); get input power for modulated and unmodulated state
   [x, y] = meshgrid(tagchar_mod.vic_pav, tagchar_mod.vic_fch);
   %     reflection coefficient unmodulated/modulated around operating point
   rho0 = interp2c(x,y, squeeze(tagchar_mod.rhoext_pav(1+tagchar_mod.settings.nf_out1+tagchar_mod.vic_fch,...
      end, tagchar_mod.vic_pav)), tagchar_mod.pos_pav, tagchar_mod.pos_fch, 'linear', NaN);
   rho1 = interp2c(x,y, squeeze(tagchar_mod.rhoext_pav(1+tagchar_mod.settings.nf_out1+tagchar_mod.vic_fch,...
        1, tagchar_mod.vic_pav)), tagchar_mod.pos_pav, tagchar_mod.pos_fch, 'linear', NaN);
   %     input power for unmodulated/modulated
   p_in0 = tagchar_mod.pav_hat * ( 1 - abs(rho0)^2 );
   p_in1 = tagchar_mod.pav_hat * ( 1 - abs(rho1)^2 );
   %     position for chip input power (unmodulated/modulated; 2Dinterp approximated by 2x1D interp)
   pos_pic0(1) = interp1(squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(1), end, :)), 1:tagchar_mod.settings.np, p_in0, 'linear');
   pos_pic0(2) = interp1(squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(2), end, :)), 1:tagchar_mod.settings.np, p_in0, 'linear');
   pos_pic1(1) = interp1(squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(1),   1, :)), 1:tagchar_mod.settings.np, p_in1, 'linear');
   pos_pic1(2) = interp1(squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(2),   1, :)), 1:tagchar_mod.settings.np, p_in1, 'linear');  
   pos_pic0 = pos_pic0(2) * mod(tagchar_mod.pos_fch,1) + pos_pic0(1) * (1-mod(tagchar_mod.pos_fch,1));
   pos_pic1 = pos_pic1(2) * mod(tagchar_mod.pos_fch,1) + pos_pic1(1) * (1-mod(tagchar_mod.pos_fch,1));
   %     chip input power for unmodulated/modulated   
   p_ic0 = interp1(1:tagchar_mod.settings.np, tagchar_mod.pic, pos_pic0, 'linear');
   p_ic1 = interp1(1:tagchar_mod.settings.np, tagchar_mod.pic, pos_pic1, 'linear');
   %     average powers taking modulation duty cycle into account (50% for EPC class-1 gen-2)
   if isempty(modsignal_fclk)
      tagchar_mod.pin_hat = p_in0;
      tagchar_mod.pic_hat = p_ic0;
   else
      % combine power, taking modulation into account
      pin_hat = p_in1 * mean(modsignal_fclk) + p_in0 * (1-mean(modsignal_fclk));
      pic_hat = p_ic1 * mean(modsignal_fclk) + p_ic0 * (1-mean(modsignal_fclk));
      % bias towards unmod to take the power supply buffer into account (until buffer is implemented)
      % 0..1, 0: mod/unmod (continuous operation), 1: unmod only (short term)
      % (for pin and pic to avoid pic > pin)
      tagchar_mod.pin_hat = p_in0 * settings.pfnc_bias + pin_hat * (1-settings.pfnc_bias);
      tagchar_mod.pic_hat = p_ic0 * settings.pfnc_bias + pic_hat * (1-settings.pfnc_bias);
   end
else
   tagchar_mod.pin_hat = NaN; % don't know
   tagchar_mod.pic_hat = NaN;
end
%     output
msg('Input power / Chip input power: Pin = %.2e uW (%.2f dBm) / Pic = %.2e uW (%.2f dBm)',...
   tagchar_mod.pin_hat*1e6, 10*log10(tagchar_mod.pin_hat)+30, tagchar_mod.pic_hat*1e6, 10*log10(tagchar_mod.pic_hat)+30);
%     variable cleanup
clear('rho0','rho1','p_in0','p_in1','pos_pic0','pos_pic1','p_ic0','p_ic1','pin_hat','pic_hat');

% is the tag functional (Pic > threshold for average mod resistance) ?
% ... this check is NOT identical to tagchar_mod.support.nonext; the nonextrapolated char. contains
%     values below pmin!
% ... for revisions below alpha 1.0 (tagchar_modulation.m below alpha 2.0)
%     tagchar_mod.support.nonext was identical to tagchar_mod.support.functional. Newer versions
%     support reflection of nonfunctional tags (char. includes power levels below operation threshold).
tagchar_mod.ind_fc = interp1(tagchar_mod.fch, 1:tagchar_mod.settings.nfch, settings.fc, 'nearest');
tagchar_mod.support.functional = tagchar_mod.support.ext && tagchar_mod.support.nonext &&...
   ~isnan(tagchar_mod.pic_hat) && tagchar_mod.picopmin(tagchar_mod.ind_fc) <= tagchar_mod.pic_hat;

% calculate power supply voltage (Vdda)
tagchar_pwr = loadmat(settings.pwrcharfile);
if ~isfield(tagchar_pwr, lower(settings.pwrchar))
   err('Unsupported power supply characteristic "%s".', settings.pwrchar);
else
   vdda = interp1(tagchar_pwr.pic, tagchar_pwr.(lower(settings.pwrchar)), 10*log10(tagchar_mod.pic_hat)+30, 'nearest', 'extrap');
end
%     output
msg('Power supply voltage: Vdda = %.2f V (saturated to bounds)', vdda);

% return some internal values
%     power levels
internal.pav = tagchar_mod.pav_hat;
internal.pin = tagchar_mod.pin_hat;
internal.pic = tagchar_mod.pic_hat;
%     functional
internal.func = tagchar_mod.support.functional;
%     power supply voltage
internal.vdda = vdda;

% decide what to do next
%     absolutely no characteristic support => modulation/reflection will not work
if ~tagchar_mod.support.ext
   msg('No characteristic support for this combination Pav,fc. Returning a zero vector.');
   modcarrier = zeros(length(carrier) - settings.grpdel_s, 1);
   return
end
%     tag is not functional
if ~tagchar_mod.support.functional
   msg('Tag is not functional for this combination incident power / carrier frequency.');
   % forced modulation/reflection (we have at least extrapolated support)
   if settings.force
      critwarn('Forced modulation/reflection.');      
   % modulation: switch to unmodulated reflection
   elseif strcmpi(settings.mode, 'modulate') && ~settings.force
      msg('Switching from modulated to unmodulated reflection.');
      modsignal_fclk = modsignal_fclk * 0;
   % filter (transmission): tag not functional and not forced => no point in calculating transmitted signal
   elseif strcmpi(settings.mode, 'filter') && ~settings.force
      modcarrier = zeros(length(carrier) - settings.grpdel_s, 1);
      return
   % to be on the safe side
   else
      criterr('Ended up in an undefined state (not functional, mode=%s, force=%i).', settings.mode, settings.force);
   end
end
%     only extrapolated support (this is only important if we continued)
if ~tagchar_mod.support.nonext
   critwarn('Only extrapolated characteristic support for this combination Pav,fc; error may be arbitrary.');
end

% output support
%     determine support boundaries
%     ... identical for all rmod because pav is calulated without modulation  => implicit assumption that
%         chip power supply is stable during modulation periods and those periods are kept short)
%     ... start/end with "no support" (1-0.1=0.9); bias towards non-NaN (-0.1); 
tagchar_mod.support.range = find(~isnan(sum(tagchar_mod.rho_pav(:, end, tagchar_mod.vic_pav), 3)));
tagchar_mod.support.bound = findzeros([0.9; isnan(sum(tagchar_mod.rho_pav(:, end, tagchar_mod.vic_pav), 3)) - 0.1; 0.9]) - 1;
%     safety check: length tagchar_mod.support.bound has to be 2,4,6,... (...no-yes-no-yes-no...)
if mod(length(tagchar_mod.support.bound), 2) ~= 0
   criterr('Characteristic support identification failed.');
end
%     output support
msg('Characteristic support:')
msg('   Extrapolated: %.1f to %.1f MHz', tagchar_mod.f(1)/1e6, tagchar_mod.f(end)/1e6)
if ~isempty(tagchar_mod.support.bound)
   mtext = sprintf('%.1f to %.1f MHz and ',...
      tagchar_mod.fch(tagchar_mod.support.bound(1:end/2))/1e6, tagchar_mod.fch(tagchar_mod.support.bound(end/2+1:end))/1e6); 
   msg('   Non-Extrapolated: %s', mtext(1:end-4));
end


% *******************************************************************************************************
% filter to obtain backscattered / transmitted signal

% check if provided channel filter coefficients are applicable, recalculate if necessary
if (settings.nfft == tagchar_mod_filterbank.settings.nfft) && (settings.nfilt == tagchar_mod_filterbank.settings.nfilt)
   coeff = tagchar_mod_filterbank.coeff;
else 
   msg('Provided filterbank channel filter coefficients not applicable; recalculating.');
   coeff = npr_coeff(settings.nfft, settings.nfilt);
end

% warning if settings.fs exceeds sampling frequency of filterbank (extrapolation has to be done)
if settings.fs * (settings.nfft-1)/(2*settings.nfft) > tagchar_mod.f(end)
   critwarn('Sampling frequency exceeds fs of characteristic, extrapolating with rho_max=%f.',...
      tagchar_mod.settings.rhomax);
end

% create overall filter impulse response w. homogenous resolution for filterbank
%     homogenous frequency resolution, unmod and mod
%     ... spline interpolation is not a good idea here, because it leads to |rho|>1 if the 
%         characteristic has been truncated to |rho|~1 (high steepness plus saturation)
[x , y ] = meshgrid(tagchar_mod.pav(tagchar_mod.vic_pav), tagchar_mod.f);
[xi, yi] = meshgrid(tagchar_mod.pav_hat, linspace(0,settings.fs/2,settings.nfft/2+1));
hh0 = interp2c(x,y, squeeze(tagchar_mod.rhoext_pav(:, end, tagchar_mod.vic_pav)), xi, yi,...
	'linear', tagchar_mod.settings.rhomax);
hh1 = interp2c(x,y, squeeze(tagchar_mod.rhoext_pav(:,   1, tagchar_mod.vic_pav)), xi, yi,...
	'linear', tagchar_mod.settings.rhomax);
%     check
if any(abs(hh0) > 1) || any(abs(hh1) > 1)
   critwarn('Resampling of characteristic failed (|rho|>1): saturating. Check sampling rate and interpolation.');
   ind0 = find(abs(hh0) > 1);
   ind1 = find(abs(hh1) > 1);
   hh0(ind0) = hh0(ind0) / abs(hh0(ind0));
   hh1(ind1) = hh1(ind1) / abs(hh1(ind1));
end

% change to transmit signal for mode filter
if strcmpi(settings.mode, 'filter')
   warn('Phase not supported in filter mode (phase of returned signal identical to carrier phase).')
   hh0 = sqrt(1 - abs(hh0).^2);
   hh1 = hh0; % just to make sure (should be unmod anyway)
end

% calculate some values
%     frequency resolution of filterbank
internal.fres = settings.fs / settings.nfft; % Hz
%     maximum error caused by filterbank not including interpolation [abs, angle]
%     ... 1/fbdiff*maximum difference between bins (cf. frequency sweep as created by plots_tag_modulation.m)
fch_range_i = [1+tagchar_mod.settings.nf_out1, tagchar_mod.settings.nf_out1 + tagchar_mod.settings.nfch];
internal.maxerr_mod = [...
   max(abs(diff(  abs(hh1(fch_range_i(1):fch_range_i(2))))))*internalsettings.fbdiff,...
   max(abs(diff(angle(hh1(fch_range_i(1):fch_range_i(2))))))*internalsettings.fbdiff*180/pi];
internal.maxerr_umd = [...
   max(abs(diff(  abs(hh0(fch_range_i(1):fch_range_i(2))))))*internalsettings.fbdiff,...
   max(abs(diff(angle(hh0(fch_range_i(1):fch_range_i(2))))))*internalsettings.fbdiff*180/pi];
%     linear model
internal.rho  = (hh1 + hh0) / 2;
internal.drho = internal.rho - hh0;

% output
%     frequency resolution
msg('Resolution of filterbank: %.3f MHz / %.3f us', internal.fres/1e6, 1e6/internal.fres);
%     and maximum error 
msg('Approximate max error for frequ between bins (highres area of tagchar, not including interpolation):')
msg('    mod: %.3e  %.3e deg', internal.maxerr_mod);
msg('   umod: %.3e  %.3e deg', internal.maxerr_umd);

% include negative frequencies (-> for npr_... functions)
% ... Note that the original npr_analysis used the IFFT, not the FFT (=> conj(hh0), conj(hh1)). The
%     modified version used here uses an FFT for the analysis (=> no modifications to hh0 and hh1 
%     necessary).
hh0 = [hh0(1:end); conj(hh0(end-1:-1:2))];
hh1 = [hh1(1:end); conj(hh1(end-1:-1:2))];

% analysis filterbank
spectrum = npr_analysis(coeff, carrier);

% free some memory (carrier needs a lot of RAM + filter will take some time)
clear('carrier');

% add initial delay to clock and compensate for modulation signal delay
clock = clock + settings.t0_s + settings.moddel_s;

% filter
%     generate index list (spectra are 50% overlapped; clock has jitter)
ind_clk = 1 + floor(interp1(clock, [1:1:length(clock)], [1:size(spectrum, 2)]*settings.nfft/2, 'linear', 'extrap'));
%     filter loop
for i = 1 : size(spectrum, 2)  
   % for t < 0 || t > modulation length: not modulated
   if (ind_clk(i) < 1) || (ind_clk(i) > length(modsignal_fclk)) || (modsignal_fclk(ind_clk(i)) == 0)
      spectrum(:, i) = spectrum(:, i) .* hh0; % not modulated
   else
      spectrum(:, i) = spectrum(:, i) .* hh1; % modulated
   end
end

% synthesis filterbank
modcarrier = npr_synthesis(coeff, spectrum);

% check imaginary part (has to be close to zero)
% ... imaginary part may be rather high if hh0 approx conj(hh1)
if mean(imag(modcarrier(settings.grpdel_s+1:end)).^2) / mean(real(modcarrier(settings.grpdel_s+1:end)).^2) > internalsettings.imagtol
   critwarn('Backscatter signal has considerable imaginary part (avg im/re = %.2g). Check filter operation for symmetry!', mean(imag(modcarrier).^2)/mean(real(modcarrier).^2));
end

% make real and remove transient (group delay); make modcarrier a column vector
modcarrier = real(modcarrier(settings.grpdel_s+1:end)).';

% multiply by abs(input impedance @ operating point) to obtain demodulator input signal in 'filter' mode
% ... reflected signal can be seen as [Vrms]@R=1Ohm (effective power: P=var(modcarrier)=U^2/R)
% ... use Zat (and not alternative) to avoid accumulation of errors (we're using mostly estimates!)
%     alternative: remove 1/Zat and add Pin/Pic
if strcmpi(settings.mode, 'filter')
   % get impedances at operating point
   %     chip input impedance
   [x , y ] = meshgrid(tagchar_mod.pic, tagchar_mod.fch);
   [xi, yi] = meshgrid(tagchar_mod.pic_hat, settings.fc);
   internal.zic = interp2c(x,y, tagchar_mod.zic, xi, yi, 'linear', NaN);
   %     assembly impedance
   internal.zat = complex(tagchar_mod.settings.rat, -1/(2*pi*settings.fc*tagchar_mod.settings.cat));
   % multiply => peak voltage on chip input impedance (no phase!)
   modcarrier = modcarrier / sqrt(real( 1/internal.zic + 1/tagchar_mod.zmod(end) + 1/internal.zat ));
end


return
% *******************************************************************************************************
% DEBUG

% % remove negative frequencies from filter frequ. response
% hh0 = [hh0(1:end/2+1)];
% hh1 = [hh1(1:end/2+1)];
% 
% figure; plot(modcarrier)
% figure; hold on; plot(carrier, 'b'); plot(modcarrier, 'r'); hold off;
% figure; hold on; plot(carrier*0.2138, 'b'); plot(modcarrier, 'r'); hold off;
%
% figure; plot(linspace(0,settings.fs/2,settings.nfft/2+1), abs(hh0))
% 
% internalsettings.demod.order     =     4; % filter order (/2 because filtfilt is used)
% internalsettings.demod.fcut      = 1/100; % normalized cutoff frequency (1: fs/2)
% internalsettings.demod.att       =    40; % dB stopband attenuation  (/2 because filtfilt is used)
% [b,a] = cheby2(internalsettings.demod.order, internalsettings.demod.att, internalsettings.demod.fcut);
% signal = filtfilt(b,a, abs(modcarrier));
% 
% figure; hold on;
% plot(signal)


%       disp('lower/upper bounds for vic_fch(1) and vic_fch(2) in pin matrix');
%       p_in0
%       squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(1), end, [floor(pos_pic0(1)), ceil(pos_pic0(1))]))
%       squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(2), end, [floor(pos_pic0(2)), ceil(pos_pic0(2))]))
%       p_in1
%       squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(1),  1, [floor(pos_pic1(1)), ceil(pos_pic1(1))]))
%       squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(2),  1, [floor(pos_pic1(2)), ceil(pos_pic1(2))]))

%    disp('reverse check for pos_pin [pin, (1)->pin, (2)->pin]');
%    [p_in0, interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(1), end, :)), pos_pic0(1), 'linear'),...
%            interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(2), end, :)), pos_pic0(2), 'linear')]
%    [p_in1, interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(1),   1, :)), pos_pic1(1), 'linear'),...
%            interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(2),   1, :)),
%            pos_pic1(2), 'linear')]   


   
%    disp('reverse check II for pos_pin [pin, (1)->pin, (2)->pin]');
%    [p_in0, interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(1), end, :)), pos_pic0, 'linear'),...
%            interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(2), end, :)), pos_pic0, 'linear')]
%    [p_in1, interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(1),   1, :)), pos_pic1, 'linear'),...
%            interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.vic_fch(2),   1, :)), pos_pic1, 'linear')]  
% 
%    [interp1(1:tagchar_mod.settings.np, tagchar_mod.pic, pos_pic0+0.5, 'linear'),...
%    interp1(1:tagchar_mod.settings.np, tagchar_mod.pic, pos_pic1+0.5, 'linear')]
%    
%    squeeze(tagchar_mod.pin(tagchar_mod.vic_fch, end, [floor(pos_pic0), ceil(pos_pic0)]))
%    squeeze(tagchar_mod.pin(tagchar_mod.ind_fch, end, round(pos_pic0)))
%    
%    squeeze(tagchar_mod.pin(tagchar_mod.vic_fch,   1, [floor(pos_pic1), ceil(pos_pic1)]))
%    squeeze(tagchar_mod.pin(tagchar_mod.ind_fch,   1, round(pos_pic1)))
%    
%    tagchar_mod.pic([floor(pos_pic0), ceil(pos_pic0), ceil(pos_pic0)+1])
%    
%    interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.ind_fch, end, :)), pos_pic0, 'linear')
%    interp1(1:tagchar_mod.settings.np, squeeze(tagchar_mod.pin(tagchar_mod.ind_fch,   1, :)), pos_pic1, 'linear')
%    
%    tagchar_mod.pin(tagchar_mod.ind_fch,   1, pos_pic1)
% 
%    disp('relerr due to fch estimate')
%    (abs(log10(tagchar_mod.fch(tagchar_mod.vic_fch)))' - abs(log10(settings.fc)) ) ./ abs(log10(settings.fc))
%       
%    disp('lower/upper bound for pic out of above positions')
%    tagchar_mod.pic([floor(pos_pic0), ceil(pos_pic0)])
%    tagchar_mod.pic([floor(pos_pic1), ceil(pos_pic1)])
%    
%    disp('these should be close [p_in, p_ic * theory]')
%    [p_in0, p_ic0 * ( abs(log10(settings.fc)) * abs(log10( tagchar_mod.rmod(end) )) * abs(log10(tagchar_mod.pav_hat)) ).^(1/4)]
%    [p_in1, p_ic1 * ( abs(log10(settings.fc)) * abs(log10( tagchar_mod.rmod(  1) )) * abs(log10(tagchar_mod.pav_hat)) ).^(1/4)]
% 
%    disp('relerr')
%    (p_ic0 * ( abs(log10(settings.fc)) * abs(log10( tagchar_mod.rmod(end) )) * abs(log10(tagchar_mod.pav_hat)) ).^(1/4) - p_in0) /p_in0
%    (p_ic1 * ( abs(log10(settings.fc)) * abs(log10( tagchar_mod.rmod(  1) )) * abs(log10(tagchar_mod.pav_hat)) ).^(1/4) - p_in1) /p_in1
% 
%       squeeze(tagchar_mod.pin(tagchar_mod.vic_fch,   1, [floor(pos_pic1), ceil(pos_pic1)]))
%       
%       tagchar_mod.pic([floor(pos_pic0), ceil(pos_pic0)])
%       
%       disp('[est min bound, est max bound, est, theo]')
%       [tagchar_mod.pic([floor(pos_pic0), ceil(pos_pic0)])', p_ic0, ...
%          p_in0 / ( abs(log10(settings.fc)) * abs(log10( tagchar_mod.rmod(end) )) * abs(log10(tagchar_mod.pav_hat)) ).^(1/4)]
%       [tagchar_mod.pic([floor(pos_pic1), ceil(pos_pic1)])', p_ic1, ...
%          p_in1 / ( abs(log10(settings.fc)) * abs(log10( tagchar_mod.rmod(1) )) * abs(log10(tagchar_mod.pav_hat)) ).^(1/4)]
%       
% 
% [settings.fc, tagchar_mod.fch(tagchar_mod.ind_fch)]
% abs(log10( tagchar_mod.rmod(end) ))




% *******************************************************************************************************
% OLD AND RUSTY

