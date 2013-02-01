% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% channel - small-scale (time-invariant)
%
% References:
% [1] K. Witrisal, A New Method to Measure Parameters of Frequency-Selective Ration Channels Using Power
%     Measurements, IEEE Transactions on Communications, Vol. 59, No. 10, October 2001, pp. 1788-1800
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
% version = channel_small()
%    Just returns the version number (string).
% [delay_s, gain, statistics] = channel_small(settings)
%    Returns gain as a function of delay (in samples) using the power-delay-profile (PDP) model in [1]. 
%    The function tries to generate the an average power-delay-profile according to settings. The final
%    results are stored in STATISTICS. The function returns a perfect average power-delay-profile 
%    if SETTINGS.DET is set to true. For SETTINGS.DET==false SETTINGS.ENSEMBLES independent (random)
%    ensembles are generated. Note that the delay axis is also randomized in this case.
%
%
% ***** Interface definition *****
% function [delay_s, gain, statistics] = channel_small(settings)
%    settings    struct containing all necessary parameters
%       .on          switch small-scale model on/off
%       .det         return the average PDP (true) or random NLOS paths with the same average PDP (false)
%       .seed        random seed for non-deterministic setup (set to NaN for normal operation)
%       .ensembles   number of independent ensembles for settings.det=false
%       .maxrays     maximum number of paths (including the LOS path, i.e. maxrays=1 means LOS only)
%       .nrays_f     (optional) force number of paths (including the LOS path)
%       .maxiter     maximum iterations for getting the rms delay spread right
%       .k           K-factor w.r.t. the LOS in dB
%       .trms        RMS delay spread in s
%       .bw          one-sided bandwidth for channel (with oversampling) in Hz
%       .fres        initial frequency resolution in Hz (might be increased to meet .trms)
%       .eps_k       tolerable relative error of K-factor (linear)
%       .eps_trms    tolerable relative error of RMS delay spread (RECOMMENDED: 10*eps_trms)
%       .fs          sampling frequency in Hz
%
%   delay_s      delays for the different paths in samples (lines are ensembles for settings.det=false)
%   gain         gain for the different paths (lines are ensembles for settings.det=false)
%   statistics   struct containing estimated statistics of the generated average PDP and warning switches
%      .av_k           estimated K-factor w.r.t. LOS in dB
%      .av_trms        estimated RMS delay spread
%      .iter           number of iterations needed to meet settings.eps_k
%      .n              number of paths (including LOS path)
%      .fres           final frequency resolution of the channel
%      .bw             final bandwidth of the channel
%      .warn_res       warning: there might be problems with the resulting frequency resolution
%      .warn_bw        warning: there might be problems with the resulting bandwidth
%      .warn_ktrms     warning: K and trms might be wrong / are off
%      .warn_maxiter   warning: maximum number of iterations reached
%      .warn_maxrays   warning: maximum number of paths reached
%      .warn_rand      warning: delay resolution insufficient to randomize (will affect narrowband corr.)
%      .corr_shift     shift of NLOS-part in order to reach trms in s
%      .corr_multf     multiplication factor for #taps and tmax in order to reach eps_trms
%      .los2nlos_f     factor in tap spacing (LOS->NLOS)/(NLOS->NLOS); (NLOS is equidist., LOS-NLOS not)
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
% - equal spacing also for LOS to first NLOS plus include trms in tmax-calculation
% - performant method to randomize delay-bins independently (would solve min-diff problem)
% ? distributing of tap energy instead of rounding (important for small number of taps)
% ? sampling of spaced-frequency correlation function (does not work so far)
%   
% *******************************************************************************************************


function [delay_s, gain, statistics] = channel_small(settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   delay_s = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% bound for tmax to avoid numerical problems in center of gravity equation (at gamma*tmax~700)
internalsettings.tg_max      = 250; % upper bound for gamma*tmax ("time constants")
internalsettings.minres_tmax =   2; % samples minimum resolution before tmax is increased in iterations


% *******************************************************************************************************
% input parameter checks

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'on','det','seed','ensembles','maxrays', 'k','trms','eps_k','eps_trms','bw','fres', 'fs'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% initialize warning switches
statistics.warn_res     = false; % there might be problems with the resulting frequency resolution
statistics.warn_bw      = false; % there might be problems with the resulting bandwidth
statistics.warn_ktrms   = false; % K and trms might be wrong / are off
statistics.warn_maxiter = false; % maximum number of iterations reached
statistics.warn_maxrays = false; % maximum number of paths reached
statistics.warn_rand    = false; % delay resolution insufficient to randomize
statistics.corr_shift   =     0; % s shift of NLOS-part in order to reach trms
statistics.corr_multf   =     0; % multiplication factor for #taps and tmax in order to reach eps_trms

% return if switched off (this should be done outside this function)
if ~settings.on
   warn('This function has been called "switched off". Is this a mistake?');
   delay_s = 0;
   gain    = 1;
   statistics.av_k    = Inf;
   statistics.av_trms =   0;
   return
end


% *******************************************************************************************************
% create average Power-Delay-Profile [1] with linear spacing meeting settings as good as possible
% ... this is like sampling the spaced-frequency-correlation-function with subsequent upsampling without
%     filtering => there will be mirrors outside statistics.bw

% K-factor: [dB] -> linear
settings.k_lin = 10^(settings.k/10);

% parameters of PDP (overall power of 0 dB ... P0=1)
los   = settings.k_lin/(settings.k_lin+1); % power of LOS component
gamma = 1/settings.trms * sqrt(2*settings.k_lin+1)/(settings.k_lin+1); % decay constant of exp. distribution
nlos  = 1/(settings.k_lin+1) * gamma; % power density NLOS

% number of taps to meet given bandwidth and frequency resolution
settings.n = ceil(2*settings.bw/settings.fres) + 1;
if settings.n > settings.maxrays
   statistics.warn_maxrays = true;
   settings.n = settings.maxrays;
end
%     force number of paths?
if isfield(settings, 'nrays_f') && ~isnan(settings.nrays_f) && settings.nrays_f >= 0
   critwarn('Forced number of paths. Note that this has priority over most other settings.');
   settings.n = settings.nrays_f;
   settings.maxiter = 1;
end
% maximum delay based on min. freq. resolution and K-factor (assymptotically correct for large #taps)
% ... RMS delay spread incorporated via approximation
eps_k_trms = min(settings.eps_k, (settings.eps_trms/3)^sqrt(2));
tmax = max(1/settings.fres*(settings.n-1)/settings.n, -log(eps_k_trms/(eps_k_trms+1))/gamma);
% tmax = max(1/settings.fres*(settings.n-1)/settings.n, -log(settings.eps_k/(settings.eps_k+1))/gamma);
if gamma*tmax > internalsettings.tg_max % truncate if necessary
   statistics.warn_ktrms = true;
   tmax = internalsettings.tg_max/gamma;
end

% iteratively generate average PDP (try to get rms delay spread right)
for iter = 1 : settings.maxiter
   %     placement of integration bounds
   intb = linspace(0, tmax, settings.n)';
   %     placement of delays at centers of gravity
   %     ... Problems (mostly for small number of taps and close to critical sampling):
   %           .) rounding (distribute tap power to adjacent taps)?
   %           .) tap spacing LOS-NLOS considerably smaller than NLOS-NLOS
   delay_s = [0; round( (statistics.corr_shift + ( exp(intb(1:end-1)*gamma) .* (1+intb(2:end)*gamma) - exp(intb(2:end)*gamma) .* (1+intb(1:end-1)*gamma) ) ./...
      (gamma * ( exp(intb(1:end-1)*gamma) - exp(intb(2:end)*gamma) )) ) * settings.fs )];
   %     generate gains
   gain = sqrt([los; -nlos/gamma * ( exp(-gamma*intb(2:end))-exp(-gamma*intb(1:end-1)) )]);  
   
   % check for identical entries in delay_s (too many rays for the given timings)
   if any(diff(delay_s)==0)
      % output a warning after the iteration is complete
      statistics.warn_res = true;
      % get a list with identical delays
      id_del = delay_s(diff(delay_s)==0);
      id_del(diff(id_del)==0) = [];
      % browse through list and combine gains for identical entries
      for i = 1 : length(id_del)
         ind = find(delay_s(:)==id_del(i));
         gain(ind(1)) = sqrt(sum(gain(ind).^2)); % geometric add.
         delay_s(ind(2:end)) = [];
         gain(ind(2:end))    = [];
      end
   end
     
   % calculate statistics
   statistics.av_k    = 10*log10(gain(1)^2/sum(gain(2:end).^2)); % K-Factor in dB
   statistics.av_trms = sqrt(var(delay_s, gain.^2))/settings.fs; % RMS delay spread
         
   % check RMS delay spread, correct if necessary
   if abs(statistics.av_trms - settings.trms) / settings.trms > settings.eps_trms      
      % there is likely a sampling problem
      %  => increase number of taps and tmax (trms asymptotically correct for settings.n>>)
      
      % break if we can't change anything any more
      if iter == settings.maxiter
         statistics.warn_maxiter = true;
         break;
      end
      if statistics.warn_maxrays
         break;
      end
      
      % increase settings.n and tmax
      %  multiplication factor for tmax and n based on relative error of eps_trms 
      corr_multf = (1+sqrt((abs(statistics.av_trms-settings.trms)/settings.trms-settings.eps_trms)/settings.eps_trms));
      %   record for statistics
      statistics.corr_multf = statistics.corr_multf * corr_multf;
      %   increase number of paths
      settings.n = ceil(settings.n * corr_multf);
      if settings.n > settings.maxrays
         statistics.warn_maxrays = true;
         settings.n = settings.maxrays;
      end
      %   increase tmax (only if resolution is critical)
      if tmax*settings.fs / settings.n < internalsettings.minres_tmax
         tmax = tmax * corr_multf;
         if gamma*tmax > internalsettings.tg_max
            tmax = internalsettings.tg_max / gamma;
         end
      end
      
   else
      break;
   end
end

% record some statistics
statistics.iter       = iter;
statistics.n          = length(delay_s);
statistics.fres       = 1 / tmax * (statistics.n-1)/statistics.n;
statistics.bw         = statistics.n * statistics.fres / 2;
statistics.los2nlos_f = diff(delay_s(1:2)) / mean(diff(delay_s(2:end)));
%     warnings
if statistics.bw < settings.bw
   statistics.warn_bw = true;
end
if statistics.fres > settings.fres
   statistics.warn_res = true;
end
if abs(statistics.av_trms - settings.trms) / settings.trms > settings.eps_trms
   statistics.warn_ktrms = true; % only necessary for trms
end


% *******************************************************************************************************
% randomize PDP (if not "deterministic")
% ... replace NLOS gains by zero-mean Gaussian RV with current variance
% ... also randomize delay_s to create random phases

if ~settings.det
   % maximum possible shift in delay
   max_diff = floor( min(diff(delay_s)) / 2 );
   if max_diff == 0 % we have a sampling problem => can't randomize delay
      statistics.warn_rand = true;
   end
   % given seed
   if ~isnan(settings.seed)
      warn('Setting RandStream (mt19937ar) seed to %s (will be reset to prior state before this fcn returns).',...
         strtrim(sprintf('%15.0f', settings.seed)));
      settings.oldDefaultStream = RandStream.getDefaultStream;
      RandStream.setDefaultStream(RandStream('mt19937ar','seed',settings.seed));
   end
   % more than one ensemble?
   if isfield(settings, 'ensembles') && settings.ensembles > 1
      gain = [gain(1)*ones(1, settings.ensembles); ...
         randn(length(gain)-1, settings.ensembles) .* repmat(gain(2:end), 1, settings.ensembles)];
      if max_diff ~= 0
         delay_s = [delay_s(1)*ones(1, settings.ensembles); ...
            repmat(delay_s(2:end),1,settings.ensembles) + randi(max_diff, [length(delay_s)-1,settings.ensembles]) - max_diff];
      else
         delay_s = repmat(delay_s,1,settings.ensembles);
      end
   else
      gain(2:end)    = randn(length(gain)-1, 1) .* gain(2:end);
      if max_diff ~= 0
         delay_s(2:end) = delay_s(2:end) + randi(max_diff, [length(delay_s)-1,1]) - max_diff;
      end
   end
   % given seed: make sure everything from here is "random" again
   if ~isnan(settings.seed)
      RandStream.setDefaultStream(settings.oldDefaultStream);
   end
end

% output warnings
if statistics.warn_res || statistics.warn_bw || statistics.warn_ktrms ||...
      statistics.warn_maxiter || statistics.warn_maxrays || statistics.warn_rand
   critwarn('Unable to meet settings. Please check statistics for details.');
   % most important details only if messages are enabled
   text = '';
   if statistics.warn_res; text=strcat(text, ', freq.res'); end
   if statistics.warn_bw; text=strcat(text, ', BW'); end
   if statistics.warn_ktrms; text=strcat(text, ', K/trms'); end
   if statistics.warn_maxiter; text=strcat(text, ', maxiter'); end
   if statistics.warn_maxrays; text=strcat(text, ', maxrays'); end
   if statistics.warn_rand; text=strcat(text, ', rand'); end
   msg('The following problems occured: %s', text(3:end));
end

% make sure delay_s is strictly monotonically rising
if any(any(diff(delay_s)<=0))
   sum(any(diff(delay_s)<=0))
   err('Delay vector is not strictly monotonically rising. Check PDP generation.');
end


end


% *******************************************************************************************************
% *******************************************************************************************************
% TESTS / OLD-AND-RUSTY

% %     check RMS delay spread
% if (settings.trms - trms_trunc(los, nlos, gamma, tmax))/settings.trms > settings.eps_trms
%    disp(sprintf('K = %.0f dB,  trms = %.1f ns', settings.k, settings.trms*1e9))
% end

% % rho, Pi, gamma, tmax
% function t2 = trms_trunc(r, p, g, t)
% gt  = g*t; % gamma * tmax
% egt = exp(-gt); % exp(-gamma*tmax)
% t2 = sqrt(  (p*( (p+g*r^2)*(2-egt.*(2+gt.*(2+gt))) - p*(egt.*(1+gt)-1).^2 )) ./ (g^2*(p+g*r^2)^2)  );
% end

% % %          settings.n = ceil(settings.n * (1+statistics.corr_ntmax));
% %          settings.n = ceil(settings.n * (1+sqrt((abs(old.av_trms-settings.trms)/settings.trms-settings.eps_trms)/settings.eps_trms)));
% %          % ...and tmax if we haven't increased it to more than internalsettings.tg_max time constants yet
% % %          tmax_test = tmax * statistics.corr_ntmax;
% % %          if gamma*tmax_test < internalsettings.tg_max
% % %             tmax = tmax_test;
% % %          else
% % %             statistics.warn_ktrms = true;
% % %          end

% % %    trms_func = @(b,g) sqrt( (1+exp(2*b*g)-exp(b*g).*(2+b.^2*g^2))./(g.^2*(exp(b*g)-1).^2) );
% % %    delay_s = [0; round( (intb(1:end-1) + trms_func(tmax/settings.n, gamma)) * settings.fs )];

%    % random delays
%    intb = sort( tmax * rand(settings.n, settings.ensembles) );
%    % placement of delays at centers of gravity (hopefully the rounding won't be much of a problem)
%    delay_s = [zeros(1,settings.ensembles); round( ( exp(intb(1:end-1,:)*gamma) .* (1+intb(2:end,:)*gamma) - exp(intb(2:end,:)*gamma) .* (1+intb(1:end-1,:)*gamma) ) ./...
%       (gamma * ( exp(intb(1:end-1,:)*gamma) - exp(intb(2:end,:)*gamma) ) ) * settings.fs )];
%    % generate gains
%    gain = sqrt([los*ones(1,settings.ensembles); -nlos/gamma * ( exp(-gamma*intb(2:end,:))-exp(-gamma*intb(1:end-1,:)) )]);


% gain_full = zeros(delay_s(end)+1,1);
% gain_full(delay_s+1) = gain;
% 
% f = [-settings.fs/2 : settings.fs/statistics.n : (statistics.n-1)/statistics.n*settings.fs/2]* statistics.n/length(gain_full)';
% sfcf = los + nlos./complex(gamma, 2*pi*f);
% 
% f_full = [-settings.fs/2 : settings.fs/length(gain_full) : (length(gain_full)-1)/length(gain_full)*settings.fs/2]';
% f_est = [-settings.bw : 2*settings.bw/length(gain) : (length(gain)-1)/length(gain)*settings.bw]';
% 
% % f_est = linspace(-settings.bw, (length(gain)-1)/length(gain)*settings.bw, 10*settings.n);
% % gain2 = [gain;zeros(9*length(gain),1)];
% 
% length(gain_full) / length(gain)
% 
% close all;
% 
% figure; hold on;
% plot(f, 10*log10(abs(sfcf)), 'b.');
% plot(f_full, 10*log10(fftshift(abs(fft(gain_full.^2)))), 'ro');
% hold off; grid on; xlim([min(f), max(f)]*1.5);
% setlabels('OLD IMPLEMENTATION, FULL', 'f [Hz]', 'abs [dB]');
% 
% figure; hold on;
% plot(f, angle(sfcf), 'b.');
% plot(f_full, fftshift(angle(fft(gain_full.^2))), 'ro');
% hold off; grid on; xlim([min(f), max(f)]*1.5);
% setlabels('OLD IMPLEMENTATION, FULL', 'f [Hz]', 'arg');
% 
% figure; hold on;
% plot(f, 10*log10(abs(sfcf)), 'b.');
% plot(f_est, 10*log10(fftshift(abs(fft(gain.^2)))), 'ro');
% hold off; grid on; xlim([min(f), max(f)]*1.5);
% setlabels('OLD IMPLEMENTATION, SPARSE', 'f [Hz]', 'abs [dB]');
% 
% figure; hold on;
% plot(f, angle(sfcf), 'b.');
% plot(f_est, fftshift(angle(fft(gain.^2))), 'ro');
% hold off; grid on; xlim([min(f), max(f)]*1.5);
% setlabels('OLD IMPLEMENTATION, SPARSE', 'f [Hz]', 'arg');
% 
% test = fftshift(fft(gain.^2));
% test2 = ifft(fftshift(test));
% 
% figure; plot(abs(test2)),
% 
% 
% % figure; hold on;
% % plot(10*log10(abs(ifft(fftshift(sfcf)))), 'b.');
% % plot(20*log10(gain), 'ro');
% % hold off; grid on; 
% % 
% % xlim([min(f), max(f)]*1.5);
% % setlabels('OLD IMPLEMENTATION, SPARSE', 'f [Hz]', 'arg');
% 
% f = [-settings.n/(2*tmax): 1/tmax : (settings.n-1)/(2*tmax)]';
% sfcf_1 = los + nlos./complex(gamma, 2*pi*f);
% % if mod(settings.n,2)==0 % make sure the inverse transformation is real
% %    g = ifft([sfcf(1:end); conj(sfcf(end-1:-1:2))]);
% % else
% %    g = ifft([sfcf(1:end); conj(sfcf(end:-1:2))]);
% % end
% 
% g = ifft(fftshift(sfcf_1));
% [sum(real(g).^2), sum(imag(g).^2)]
% g = abs(g);
% 
% figure; plot(g, 'b.')
% 
% d = round([1/settings.fs : 1/(2*settings.bw) : tmax]' * settings.fs) - 1;
% 
% 
% % figure; hold on;
% % plot(f, abs(sfcf), 'b.');
% % plot(f_full, fftshift(abs(fft(gain_full.^2))), 'ro');
% % plot([settings.bw, settings.bw], [min(abs(sfcf)), max(abs(sfcf))], 'g-');
% % hold off; grid on; xlim([min(f), max(f)]*1.5);
% % setlabels('OLD IMPLEMENTATION, FULL', 'f [Hz]', 'abs');
% 
% 
% g_full = zeros(d(end)+1,1);
% g_full(d+1) = g;
% 
% f_est = [-settings.fs/2 : settings.fs/length(g_full) : (length(g_full)-1)/length(g_full)*settings.fs/2]';
% 
% figure; hold on;
% plot(f , abs(sfcf), 'b.');
% plot(f_est, abs(fft(g_full.^2)), 'ro');
% hold off; grid on; xlim([min(f), max(f)]);
% 
% figure; hold on;
% plot(f , angle(sfcf), 'b.');
% plot(f_est, angle(fft(g_full.^2)), 'ro');
% hold off; grid on; xlim([min(f), max(f)]);
% 
% 
% [settings.k, 10*log10(g(1)^2/sum(g(2:end).^2))]
% [settings.trms, sqrt(var(d/settings.fs, g.^2))]
% 
% close all; 
% 
% 
% figure; hold on;
% plot(d, 20*log10(g), 'b.');
% plot(delay_s, 20*log10(gain), 'rx');
% hold off; grid on;
% 
% % n = [1 : 1 : delay_s(end)] / settings.fs;
% % g = ( complex(0, nlos*expint(n*gamma)) - expint(n*complex(gamma,2*pi)) ) / (2*pi);
% % g = sqrt([los, abs(g)]);

