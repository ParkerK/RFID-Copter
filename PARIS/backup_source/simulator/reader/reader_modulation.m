% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% reader - modulation (plus PIE)
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
% version = reader_modulation()
%    Just returns the version number (string).
% settings  = reader_modulation(data, settings)
%    Returns completed and sanitized settings (e.g. the needed carrier length to encode and modulate
%    DATA according to SETTINGS).
% carrier = reader_modulation(carrier, data, settings)
%    Returns the carrier modulated by DATA according to SETTINGS and starting at time SETTINGS.T0. Uses 
%    cosine rolloff for band limitation. No additional filtering is done. Checks SETTINGS for compliance  
%    with standard and corrects them if necessary (warning issued in case of corrections). 
%
%
% ***** Interface definition *****
% function carrier = reader_modulation(carrier, data, settings)
%    carrier    (optional) carrier signal 
%    data       data bits to modulate (alphabet {0,1} required!)        
%    settings   struct containing all necessary parameters
%       .forcesettings   if true, no consistency/conformity checks are done (violation of norm possible)
%       .modulation      modulation type {'DSB-ASK', 'SSB-ASK', 'PR-ASK'}
%       .leadin          {'preamble', 'frame-sync' or 'framesync'}
%       .tari            duration of data-0 in s
%       .lf              backscatter link frequency in Hz {40...640 kHz}
%       .dr              divide ratio DR {64/3, 8}
%       .delimiter       preamble/frame-sync initial delimiter duration in s ; -1: random
%       .x               additional duration of data-1 compared to data-0 (times tari) ; -1: random
%       .pw              RF pulsewidth (times tari) ; -1: random
%       .moddepth        modulation depth in percent ; - 1: random
%       .trf             rise=fall (times tari) ; -1: random
%       .t0              s of unmodulated carrier before and after modulation
%       .fs              sampling (HF) frequency in Hz
%
%    carrier     modulated carrier or returned (sanitized, completed) settings
%    -> additional settings     
%       .rtcal         preamble/framesync: length of data-0 plus data-1 in s
%       .trcal         preamble: specifies backscatter link frequency in s
%       .tari_s        .tari in samples
%       .x_s           .x in samples
%       .pw_s          .pw in samples
%       .trf_s         .trf in samples
%       .t0_s          .t0 in samples
%       .trcal_s       .trcal in samples
%       .rtcal_s       .rtcal in samples
%       .delimiter_s   .delimiter in samples
%       .leadin_s      length of leadin (preamble or framesync => see .preamble) in samples
%       .length_s      length of needed carrier signal to modulate in samples
%       .length        length of needed carrier signal to modulate in s
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
% = fix overshoots for SSB-ASK
% ? slight changes in moddepth, PW and x during transmission
% ? PR-ASK also for data vectors with odd length required (problem with carrier phase after modulation)
%
% *******************************************************************************************************

function carrier = reader_modulation(varargin)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   carrier = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings
internalsettings.tari_min      =   6.25e-6; % s
internalsettings.tari_max      =     25e-6; % s
internalsettings.trcal_min     =       1.1; % times rtcal
internalsettings.trcal_max     =         3; % times rtcal
internalsettings.delimiter_min = 11.875e-6; % s (or pw_min ... will be calculated later)
internalsettings.delimiter_max = 13.125e-6; % s
internalsettings.x_min         =       0.5; % times tari
internalsettings.x_max         =         1; % times tari
internalsettings.moddepth_min  =        80; % percent
internalsettings.moddepth_typ  =        90; % percent
internalsettings.moddepth_max  =       100; % percent
internalsettings.pw_min1       =     0.265; % times tari
internalsettings.pw_min2       =      2e-6; % s
internalsettings.pw_max        =     0.525; % times tari
internalsettings.trf_min       =         0; % times tari 
internalsettings.trf_max       =      0.33; % times tari (just dummy, will be set to PW/2)
internalsettings.t0_min        =         0; % s
internalsettings.t0_max        =       Inf; % s
internalsettings.trcal_tol     =      1e-6; % times tari tolerance for trcal
internalsettings.lf_tol        =      1e-3; % Hz tolerance for LF


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% number of input parameters
if nargin < 2
   criterr('Not enough input arguments.');
elseif nargin == 2
   justreturnsettings = true;
   data     = varargin{1}(:);
   settings = varargin{2};
elseif nargin == 3
   justreturnsettings = false;
   carrier  = varargin{1}(:);
   data     = varargin{2}(:);
   settings = varargin{3};
elseif nargin > 3
   criterr('Too many input arguments.');
end

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'forcesettings', 'modulation', 'leadin', 'tari', 'lf', 'dr', 'delimiter',...
   'x', 'pw', 'moddepth', 'trf', 't0', 'fs'};
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

% check if data length is even when PR-ASK is selected
if strcmpi(settings.modulation, 'pr-ask') && mod(length(data),2)~=0
   err('Length of data has to be even if PR-ASK is selected.')
end

% check DR
if settings.dr ~= 8 && settings.dr ~= 64/3
   err('Only divide ratios DR = 64/3 or 8 allowed.');
end


% *******************************************************************************************************
% check (an if necessary sanitize) settings

if settings.forcesettings
   critwarn('Forced setup. No consistency/conformity checks are done.');
   settings.rtcal = (2 + settings.x) * settings.tari; % s
   settings.trcal = settings.dr/settings.lf; % s
else
   % Tari (length of data-0 symbol)
   settings.tari = sanitize_setting('Tari', settings.tari, [internalsettings.tari_min, internalsettings.tari_max]);
   
   % x (length of data-1 minus length of data-0)
   if settings.x >= 0 % manual mode ... correct if necessary
      settings.x = sanitize_setting('x', settings.x, [internalsettings.x_min, internalsettings.x_max]);
   else % automatic mode (uniform distribution between min and max)
      settings.x = internalsettings.x_min + rand * (internalsettings.x_max-internalsettings.x_min);
   end

   % modulation depth 
   if settings.moddepth >= 0 % manual mode ... correct if necessary
      settings.moddepth = sanitize_setting('modulation depth', settings.moddepth, [internalsettings.moddepth_min, internalsettings.moddepth_max]);
   else % automatic mode (Gaussian around typ with (min-max)/2 = 3 sigma, truncated to min and max)
      % prevent outliers created by Gaussian distribution
      % (will not be executed again if first execution settings.moddepth = ... is within borders -> like do-while)
      while (settings.moddepth < internalsettings.moddepth_min) || (settings.moddepth > internalsettings.moddepth_max)
         settings.moddepth = internalsettings.moddepth_typ + randn * sqrt((internalsettings.moddepth_max-internalsettings.moddepth_min)/6);
      end
      % just to make sure (changes in behavior of while statement?)
      if (settings.moddepth < internalsettings.moddepth_min) || (settings.moddepth > internalsettings.moddepth_max)
         criterr('Automatic (random) setup of modulation depth failed. This means the behavior of while-statements has been changed (Matlab).')
      end
   end
   
   % pulse width PW
   %     Tari is ok now => calculate pw_min
   internalsettings.pw_min = max(internalsettings.pw_min1, internalsettings.pw_min2/settings.tari);
   if settings.pw >= 0 % manual mode ... correct if necessary
      settings.pw = sanitize_setting('pulsewidth PW', settings.pw, [internalsettings.pw_min, internalsettings.pw_max]);
   else % automatic mode (uniform distribution between min and max)
      settings.pw = internalsettings.pw_min + rand * (internalsettings.pw_max-internalsettings.pw_min);
   end
   
   % length of delimiter (preamble and frame-sync)
   if settings.delimiter >= 0 % manual mode ... correct if necessary
      settings.delimiter = sanitize_setting('delimiter', settings.delimiter, [internalsettings.delimiter_min, internalsettings.delimiter_max]);
   else % automatic mode (uniform distribution between min and max)
      settings.delimiter = internalsettings.delimiter_min + rand*(internalsettings.delimiter_max-internalsettings.delimiter_min);
   end

   % rise/fall time
   %     PW and delimiter are ok now => calculate tfr_max (maximum PW or (1-PW) or delimiter width)
   internalsettings.trf_max = min([internalsettings.trf_max, settings.pw, (1-settings.pw), settings.delimiter/settings.tari]);
   if settings.trf >= 0 % manual mode ... correct if necessary
      settings.trf = sanitize_setting('rise/fall time', settings.trf, [internalsettings.trf_min, internalsettings.trf_max]);
   else % automatic mode (uniform distribution between min and max)
      settings.trf = internalsettings.trf_min + rand*(internalsettings.trf_max-internalsettings.trf_min);
   end 
   
   % calculate RTcal (Tari and X now ok)
   settings.rtcal = (2 + settings.x) * settings.tari; % s
   
   % calculate TRcal out of LF and DR
   %     sanitize with tolerance to avoid unnecessary warnings in case this fcn is run twice 
   %     (once to get the necessary carrier length and once for modulation)
   settings.trcal = sanitize_setting('TRcal', settings.dr/settings.lf,...
      [internalsettings.trcal_min, internalsettings.trcal_max]*settings.rtcal, internalsettings.trcal_tol*settings.tari);
   
   % calculate LF with new TRcal
   %     sanitize with tolerance to avoid unnecessary warnings in case this fcn is run twice 
   %     (once to get the necessary carrier length and once for modulation)
   lf = settings.dr / settings.trcal;
   settings.lf = sanitize_setting('LF', settings.lf, [lf, lf], internalsettings.lf_tol);
   
   % time delay t0
   settings.t0 = sanitize_setting('t0', settings.t0, [internalsettings.t0_min, internalsettings.t0_max]);
end

% if PR-ASK is selected, moddepth is always 100%
% => set to 100%, otherwise comparison of settings to measurements will fail
if strcmpi(settings.modulation, 'pr-ask')
   settings.moddepth = 100;
end


% *******************************************************************************************************
% discretize time (possible violation of ISO neglected because .tari, .x, .pw >> 1/fs)

settings.tari_s      = round(settings.tari * settings.fs);
settings.x_s         = round(settings.x * settings.tari  * settings.fs);
settings.pw_s        = round(settings.pw * settings.tari * settings.fs);
settings.trf_s       = round(settings.trf * settings.tari * settings.fs);
settings.t0_s        = round(settings.t0 * settings.fs);
settings.trcal_s     = round(settings.trcal * settings.fs);
settings.rtcal_s     = round(settings.rtcal * settings.fs);
settings.delimiter_s = round(settings.delimiter * settings.fs);

% make sure trf_s is even (because will be divided by 2 below)
if mod(settings.trf_s,2) ~= 0
   settings.trf_s = settings.trf_s - 1;
end

% length of preamble or frame-sync => lead-in
switch lower(settings.leadin)
   case 'preamble'
      settings.leadin_s = settings.trf_s/2 + settings.delimiter_s + settings.tari_s + settings.rtcal_s + settings.trcal_s;
   case {'frame-sync', 'framesync'}
      settings.leadin_s = settings.trf_s/2 + settings.delimiter_s + settings.tari_s + settings.rtcal_s;
   otherwise
      err('Unsupported lead-in "%s%', settings.leadin);
end

% define length of modulated part (whole bitstream: complete) at sampling frequency
settings.length_s = settings.leadin_s +...                         % preamble, frame-sync (or nothing at all)
                    (length(data)-sum(data))*settings.tari_s +...  % data-0 symbols
                    sum(data)*(settings.tari_s+settings.x_s) +...  % data-1 symbols
                    settings.trf_s/2;                              % "fade-out"

% just return length if requested
if justreturnsettings
   settings.length_s = settings.length_s + 2*settings.t0_s; % add leading/trailing t0 [samples]
   settings.length   = settings.length_s / settings.fs; % [s]
   carrier = settings;
   return
end
               
% check length of carrier signal (now that timing is fixed)
if length(carrier) < 2*settings.t0_s + settings.length_s
   err('Insufficient length of carrier signal for given data vector (%i < %i).', length(carrier), 2*settings.t0_s + settings.length_s);
end


% *******************************************************************************************************
% create modulation signal
%     The fact that trf is not 90%-10%, but 100-0% here should allow for some filtering later-on
%     without violating the norm.

% define amplidute values (min/max modulation values and phase reversal multiplier)
switch upper(settings.modulation)
   case {'DSB-ASK', 'SSB-ASK'}
      am_0 = 1-settings.moddepth/100;
      am_1 = 1;
      % phase change after each bit (none for DSB/SSB-ASK)consignment
      phase_reversal_multiplier = 1;
   case 'PR-ASK'
      % min/max modulation values
      am_1 = 1; % plus/minus
      am_0 = 0;
      % phase change after each bit (only for PR-ASK)
      phase_reversal_multiplier = -1;
   otherwise
      err('Unsupported modulation type: %s (only SSB-ASK, DSB-ASK or PR-ASK allowed).', settings.modulation)
end

% define rolloffs etc for different modulation schemes
if settings.trf_s > 0 % settings.trf_s is always even (check above)
   switch upper(settings.modulation)
      case {'DSB-ASK', 'SSB-ASK'}
         % windowed rolloffs
         rolloff_start = (sin(linspace(     0,   pi/2, settings.trf_s/2))'+1)/2 * (am_1-am_0) + am_0;
         rolloff_mid   = (sin(linspace(  pi/2, 3*pi/2, settings.trf_s  ))'+1)/2 * (am_1-am_0) + am_0;
         rolloff_end   = (sin(linspace(3*pi/2, 2*pi  , settings.trf_s/2))'+1)/2 * (am_1-am_0) + am_0;
         % rolloffs for end of bitstream (for smooth transition)
         rolloff_last  = (sin(linspace(3*pi/2, 5*pi/2, settings.trf_s))'+1)/2*(am_1-am_0) + am_0;
      case 'PR-ASK'
         % windowed rolloffs
         rolloff_start = (sin(linspace(   0,   pi/2, settings.trf_s/2))'+1)/2 * am_1;
         rolloff_mid   = (sin(linspace(pi/2, 3*pi/2, settings.trf_s  ))'+1)/2 * am_1;
         rolloff_end   = (sin(linspace(pi/2,     pi, settings.trf_s/2))'-1)/2 * am_1;
         % rolloffs for start/end of bitstream (for smooth transition)
         rolloff_last  = (sin(linspace(3*pi/2, 5*pi/2, settings.trf_s))+1)/2' * am_1;
   end
else
   rolloff_start = [];
   rolloff_mid   = [];
   rolloff_end   = [];
   rolloff_last  = [];
end

% data-0 modulation signal
data_0 = [rolloff_start;...
   ones(settings.tari_s - settings.pw_s - settings.trf_s, 1) * am_1;...
   rolloff_mid;...
   ones(settings.pw_s - settings.trf_s, 1) * am_0;...
   rolloff_end];

% data-1 modulation signal
data_1 = [rolloff_start;...
   ones(settings.tari_s + settings.x_s - settings.pw_s - settings.trf_s, 1) * am_1;...
   rolloff_mid;...
   ones(settings.pw_s - settings.trf_s, 1) * am_0;...
   rolloff_end];

% preamble, frame-sync or nothing => lead-in
   switch upper(settings.modulation)
      
      case {'DSB-ASK', 'SSB-ASK'}
         switch lower(settings.leadin)
            case 'preamble'
               leadin = [rolloff_mid;...                                                % delimiter
                  ones(settings.delimiter_s - settings.trf_s, 1) * am_0;...             %   :
                  rolloff_end;...                                                       % delimiter
                  data_0;...                                                            % data-0
                  rolloff_start;...                                                     % RTcal
                  ones(settings.rtcal_s - settings.pw_s - settings.trf_s, 1) * am_1;... %   :
                  rolloff_mid;...                                                       %   :
                  ones(settings.pw_s - settings.trf_s, 1) * am_0;...                    %   :
                  rolloff_end;                                                          % RTcal
                  rolloff_start;...                                                     % TRcal
                  ones(settings.trcal_s - settings.pw_s - settings.trf_s, 1) * am_1;... %   :
                  rolloff_mid;...                                                       %   :
                  ones(settings.pw_s - settings.trf_s, 1) * am_0;...                    %   :
                  rolloff_end];                                                         % TRcal
            case {'frame-sync', 'framesync'}
               leadin = [rolloff_mid;...                                                % delimiter
                  ones(settings.delimiter_s - settings.trf_s, 1) * am_0;...         %   :
                  rolloff_end;...                                                       % delimiter
                  data_0;...                                                            % data-0
                  rolloff_start;...                                                     % RTcal
                  ones(settings.rtcal_s - settings.pw_s - settings.trf_s, 1) * am_1;... %   :
                  rolloff_mid;...                                                       %   :
                  ones(settings.pw_s - settings.trf_s, 1) * am_0;...                    %   :
                  rolloff_end];                                                         % RTcal
         end
         phase_reversal = 1; % initial phase zero deg for DSB/SSB-ASK
         
      case 'PR-ASK'
         switch lower(settings.leadin)
            case 'preamble'
               leadin = [rolloff_mid;...                                                % delimiter
                  ones(settings.delimiter_s - settings.trf_s, 1) * am_0;...             %   :
                  rolloff_end;...                                                       % delimiter
                  -data_0;...                                                            % data-0
                  rolloff_start;...                                                     % RTcal
                  ones(settings.rtcal_s - settings.pw_s - settings.trf_s, 1) * am_1;... %   :
                  rolloff_mid;...                                                       %   :
                  ones(settings.pw_s - settings.trf_s, 1) * am_0;...                    %   :
                  rolloff_end;                                                          % RTcal
                  -rolloff_start;...                                                     % TRcal
                  -ones(settings.trcal_s - settings.pw_s - settings.trf_s, 1) * am_1;... %   :
                  -rolloff_mid;...                                                       %   :
                  -ones(settings.pw_s - settings.trf_s, 1) * am_0;...                    %   :
                  -rolloff_end];                                                         % TRcal
               phase_reversal = 1; % start w. non-reversed phase
               
            case {'frame-sync', 'framesync'}
               leadin = [rolloff_mid;...                                                % delimiter
                  ones(settings.delimiter_s - settings.trf_s, 1) * am_0;...         %   :
                  rolloff_end;...                                                       % delimiter
                  -data_0;...                                                            % data-0
                  rolloff_start;...                                                     % RTcal
                  ones(settings.rtcal_s - settings.pw_s - settings.trf_s, 1) * am_1;... %   :
                  rolloff_mid;...                                                       %   :
                  ones(settings.pw_s - settings.trf_s, 1) * am_0;...                    %   :
                  rolloff_end];                                                         % RTcal
               phase_reversal = -1; % start w. reversed phase
         end
   end

% length of data-0 and data-1 signals (just for convenience)
length_0       = settings.tari_s;
length_1       = settings.tari_s + settings.x_s;

% prepare array for modulation data (might be very large vector)
modsignal = ones(settings.length_s, 1);

% place lead-in
modsignal(1:settings.leadin_s) = leadin;

% create signal bit per bit
%     helper variables
index = 1 + settings.leadin_s; % write index in modulation signal
%     create
for i = 1 : length(data)
   % place bit signal
   if data(i) == 0
      modsignal(index:index+length_0-1) = data_0 * phase_reversal;
      index = index + length_0;
   else
      modsignal(index:index+length_1-1) = data_1 * phase_reversal;
      index = index + length_1;
   end
   % reverse phase for PR-ASK
   phase_reversal = phase_reversal * phase_reversal_multiplier;
end

% shift index to end of signal (last rolloff lasts longer)
index = index + settings.trf_s/2;

% just a small check if this is really the predicted end
if index ~= settings.length_s + 1
   criterr('Did not reach the expected length of modulation signal.')
end

% replace rolloff at end of signal
modsignal(end-settings.trf_s+1:end) = rolloff_last;


% *******************************************************************************************************
% modulate

switch upper(settings.modulation)
   case {'DSB-ASK', 'PR-ASK'}
      % use both sidebands
      carrier(1+settings.t0_s:settings.t0_s+settings.length_s) =...
         modsignal .* carrier(1+settings.t0_s:settings.t0_s+settings.length_s);
   case 'SSB-ASK'
      critwarn('SSB-ASK is not fully supported (produces overshoots)!');
      % use only LSB
      carrier(1+settings.t0_s:settings.t0_s+settings.length_s) =...
         modsignal .* carrier(1+settings.t0_s:settings.t0_s+settings.length_s) +...
         imag(hilbert(modsignal)) .* imag(hilbert(carrier(1+settings.t0_s:settings.t0_s+settings.length_s)));
   otherwise
      criterr('Unsupported modulation. At least one other check for this was already passed!'); 
end



% % **********
% % DEBUG
% figure; hold on;
% plot([0:length(data_0)-1]/settings.fs, data_0,'b');
% plot([0:length(data_1)-1]/settings.fs, data_1, 'k'); xlabel('\mu s');
% plot([0:length(data_1)-1]/settings.fs, ones(size(data_1))*(am_1+am_0)/2 ,'g'); xlabel('s');
% hold off; grid on; 
% 
% figure; hold on;
% plot([0:length(leadin)-1]/settings.tari_s, leadin,'r');
% plot([0:length(data_0)-1]/settings.tari_s, data_0,'b');
% plot([0:length(data_1)-1]/settings.tari_s, data_1, 'k');
% plot([0:length(data_1)-1]/settings.tari_s, ones(size(data_1))*(am_1+am_0)/2,'g'); xlabel('x Tari');
% hold off; grid on;
% 
% figure;
% plot([0:length(modsignal)-1]/settings.tari_s, modsignal,'b');
% % **********
