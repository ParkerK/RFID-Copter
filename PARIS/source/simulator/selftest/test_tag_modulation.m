% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: tag - modulation
%
% ATTENTION: results_tag_modulation.m is not able to model the behavior of tag_modulation.m completely
%            (especially functionality checks) => it may happen that this function crashes with index
%            problems. In that case rerun results_tag_modulation to create a new set of random setups.
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
% version = test_tag_modulation()
%    Just returns the version number (string).
% sumoferrors = test_tag_modulation(sumoferrors)
%    Tests the function tag_modulation, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_tag_modulation(sumoferrors);
%    sumoferrors    sum of errors found not including this fcn call
%
%    sumoferrors    sum of errors found including this fcn call (+1 for each erroneous tested
%                   functionality)
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
% ? include tsub_s, tinv_s, tvio_s checks ("pulsewidth" tests)
%
% *******************************************************************************************************


function sumoferrors = test_tag_modulation(sumoferrors)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   sumoferrors = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% tolerances [avg, std-dev] or avg
%     timings
%        ... plus nfft for LF: can be nfft/2 @ start and nfft/2 @ end of subbit due to clock jitter
%        ... plus tclk for trailing t0 (hard to predict because of clock jitter); length should be <= predicted
%        ... clock (phase) shift is corrected for t0
internalsettings.tol_clk =             2; % samples tolerance for tag clock length (jitter!)
internalsettings.tol_lf  = [0.01, 0.001]; % [avg, std-dev] LF frequency (ratio of true LF)
internalsettings.tol_t0  =          0.05; % leading/trailing delay (ratio of true value)
%     functionality of tag_power (impedances, power levels, power supply voltage)
internalsettings.tol_zic =          0.01; % chip input impedance Zic
internalsettings.tol_pwr =          0.05; % power levels (pav, pin, pic) and pwr supply voltage (vdda)
%     modulation levels (magnitude and phase)
internalsettings.tol_abs = [0.01, 0.005]; % backscatter amplitude (std.dev: not much data!)
internalsettings.tol_arg = [0.10, 0.050]; % deg backscatter phase (std.dev: not much data!)

% "10%" and "90%" thresholds for linktiming_tag
internalsettings.linktiming_tag.levels = [0.4, 0.6]; % set rel. close to 0.5 to be robust against overshoots

% threshold for phase jumps before phase is rewrapped
internalsettings.pjmp_th = 0.85 * pi;


% *******************************************************************************************************
% initialization

% output
disp('   = tag_modulation **');
disp(sprintf('     checking power levels and supply voltage (pav, pin, pic, vdda) +/-%g%%',...
   internalsettings.tol_pwr*100));
disp(sprintf('     checking [avg +/-, std <=]: magnitude [%g, %g], phase [%g, %g]deg, LF [%g%% + nfft, %g%%]',...
   internalsettings.tol_abs, internalsettings.tol_arg, internalsettings.tol_lf*100));
disp(sprintf('     checking leading t0 +/-%g%% + nfft/2 and trailing t0 +/-%g%% + nfft/2 + tclk (jitter)',...
   internalsettings.tol_t0*100, internalsettings.tol_t0*100));

% load settings and expected results
data = loadmat('results_tag_modulation.mat', '     ');
settings     = data.settings;
results      = data.results;
partitioning = data.partitioning;


% *******************************************************************************************************
% run tests and evaluate

index = 1; % current index in partitioning.indices

for i = 1 : length(settings)

   % once for each new test block
   if i == partitioning.indices(index) + 1
      % output
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
      if strcmpi(partitioning.names{index}, 'filter')
         disp(sprintf('        additional check: checking chip input impedance Zic +/-%g%%', internalsettings.tol_zic*100));
      end  
      % reset errors and errortext
      errors    = 0;
      errortext = '';
      % step index
      index = index + 1;
   end
   
   %  re-initialize carrier if requested (e.g. at beginning of each block
   if settings{i}.newcarrier;
      basecarrier = reader_oscillator(settings{i}.reader.oscillator);
   end
     
   % prepare clock and carrier
   %     tag clock in case of modulation (no clock for "filter")
   if strcmpi(settings{i}.mode, 'modulate')
      % let tag_modulation set the appropriate lengths (important for trailing t0, also a tested fcn)
      set_mod = tag_modulation(settings{i}.data, settings{i}.tag.modulation);
      settings{i}.clen_s = set_mod.length_s; % carrier length
      settings{i}.tag.modulation.t0 = set_mod.t0_s / settings{i}.fs; % leading/trailing unmod. part
      % set tag clock length (samples @ fclk) and create clock
      settings{i}.tag.clock.length = set_mod.length_fclk_s;
      clock = tag_clock(settings{i}.tag.clock);      
   else
      clock = []; % no clock needed in filter mode
      settings{i}.clen_s = length(basecarrier); % use full carrier length
   end
   %     select carrier snippet and set power (plus additional lenght)
   carrier = basecarrier(1:settings{i}.clen_s + settings{i}.addclen_s) * sqrt(settings{i}.pav);  
   
   % modulate (cut carrier to clock length to speed up the process)
   [modulated, internal] = tag_modulation(carrier, settings{i}.data, clock, settings{i}.tag.modulation);
      
   % phase/amplitude estimation
   param_carrier   = est_sinusoid( carrier(1:length(modulated)), settings{i}.other.est_sinusoid);
   param_modulated = est_sinusoid(modulated,                     settings{i}.other.est_sinusoid);
   phase     = (param_modulated.p - param_carrier.p) * pi/180;
   magnitude = param_modulated.a ./ param_carrier.a;
   %     if phase contains jumps close to pi (and only then!): try to "wrap" phase
   if max(abs(diff(phase))) > internalsettings.pjmp_th
      phase = mod(phase, pi);
   end
      
   % timing detection if modulated ("if we have data")
   if sum(results{i}.data) ~= 0
      timings = linktiming_tag(modulated, internalsettings.linktiming_tag);
   else
      timings = struct('tran_s', [1, length(modulated)]);
   end
   
   % average magnitude/phase
   %     modulated data expected
   if sum(results{i}.data) ~= 0
      %     create grid with subbit positions, use estimated LF to avoid conflicts with LF
      bitpos = round( (timings.tran_s(1) + (0.5+[0:1:length(results{i}.data)-1]) * timings.tlf_s/2) ...
         / (settings{i}.other.est_sinusoid.n*settings{i}.other.est_sinusoid.ol/100) );
      %     indices of modulated and unmodulated parts
      ind_mod = find(results{i}.data == 1);
      ind_umd = find(results{i}.data == 0);
   %     no data => use whole vector
   else
      bitpos  = 1:length(magnitude);
      ind_mod = 1:length(bitpos);
      ind_umd = 1:length(bitpos);
   end
   
   % average magnitude/phase (always create both vectors)
   %     no characteristic support =>  should return a  zero vector => estimated phase will be NaN
   if ~results{i}.func.pav
      est_mod = [mean(magnitude), std(magnitude)];
      est_umd = [mean(magnitude), std(magnitude)];
   %     created results{i} should take care of all other cases
   else 
      est_mod = [mean(magnitude(bitpos(ind_mod))) * exp(complex(0, mean(phase(bitpos(ind_mod))))),...
                  std(magnitude(bitpos(ind_mod))) * exp(complex(0,  std(phase(bitpos(ind_mod)))))];
      est_umd = [mean(magnitude(bitpos(ind_umd))) * exp(complex(0, mean(phase(bitpos(ind_umd))))),...
                  std(magnitude(bitpos(ind_umd))) * exp(complex(0,  std(phase(bitpos(ind_umd)))))];
   end
   
   % align phase (remove phase ambiguities created by est_sinusoid.m)
   %     unmodulated
   if abs(angle(est_umd(1)) - mod(angle(results{i}.gain(1)), pi)) > internalsettings.pjmp_th
      est_umd(1) = -est_umd(1);
   end
   %     modulated
   if abs(angle(est_mod(1)) - mod(angle(results{i}.gain(2)), pi)) > internalsettings.pjmp_th
      est_mod(1) = -est_mod(1);
   end  
   
   % if in filter mode: divide by expected Rin to obtain filter gain factor again (before all that)
   % ... no characteristic support => estimated power levels will be NaN
   % ... Attention: For the test characteristic there is no connection between the impedances and
   %                the power levels (unlike for the real char). Hence using * sqrt(real(1/Zin)) here
   %                and Zic||rmod*Pic/Pin in tag_modulation won't work.
   if strcmpi(settings{i}.mode, 'filter') && results{i}.func.pav
      est_mod = est_mod * sqrt(real(1/results{i}.zin_filter)); % can't hurt
      est_umd = est_umd * sqrt(real(1/results{i}.zin_filter));
   end
   
   % compare to expected results
   %     check set clock length
   [errors, errortext] = check_result(i, length(clock), results{i}.clk_len,...
      errors, errortext, 'clk_len', 'relerr', internalsettings.tol_clk);
   %     check power levels and supply voltage (only if possible)
   [errors, errortext] = check_result(i, internal.pav, results{i}.pav,...
      errors, errortext, 'pav', 'relerr', internalsettings.tol_pwr);
   if results{i}.func.pav % only if we have at least extrapolated characteristic support
      [errors, errortext] = check_result(i, internal.pin, results{i}.pin,...
         errors, errortext, 'pin', 'relerr', internalsettings.tol_pwr);
      [errors, errortext] = check_result(i, internal.pic, results{i}.pic,...
         errors, errortext, 'pic', 'relerr', internalsettings.tol_pwr);
      [errors, errortext] = check_result(i, internal.vdda, results{i}.vdda,...
         errors, errortext, 'vdda', 'relerr', internalsettings.tol_pwr);
   end
   %     check magnitude/phase: unmodulated
   [errors, errortext] = check_result(i,   abs(est_umd(1)), abs(results{i}.gain(1)),...
      errors, errortext, 'm0(avg)', 'mse', internalsettings.tol_abs(1).^2);
   [errors, errortext] = check_result(i,   abs(est_umd(2)), 0,...
      errors, errortext, 'm0(std)', 'range', [-Inf, internalsettings.tol_abs(2)]);
   [errors, errortext] = check_result(i, angle(est_umd(1)), mod(angle(results{i}.gain(1)), pi),...
      errors, errortext, 'p0(avg)', 'mse', internalsettings.tol_arg(1).^2);
   [errors, errortext] = check_result(i, angle(est_umd(2)), 0,...
      errors, errortext, 'p0(std)', 'range', [-Inf, internalsettings.tol_arg(2)]);
   %     only in case of modulation
   if sum(results{i}.data) ~= 0
      %     check magnitude/phase: modulated
      [errors, errortext] = check_result(i,   abs(est_mod(1)), abs(results{i}.gain(2)),...
         errors, errortext, 'm1(avg)', 'mse', internalsettings.tol_abs(1).^2);
      [errors, errortext] = check_result(i,   abs(est_mod(2)), 0,...
         errors, errortext, 'm1(std)', 'range', [-Inf, internalsettings.tol_abs(2)]);
      [errors, errortext] = check_result(i, angle(est_mod(1)), mod(angle(results{i}.gain(2)), pi),...
         errors, errortext, 'p1(avg)', 'mse', internalsettings.tol_arg(1).^2);
      [errors, errortext] = check_result(i, angle(est_mod(2)), 0,...
         errors, errortext, 'p1(std)', 'range', [-Inf, internalsettings.tol_arg(2)]);
      %     check timings (no relerr here because of additional nfft, etc.)
      [errors, errortext] = check_result(i, timings.tlf_s, settings{i}.fs/results{i}.lf,...
         errors, errortext, 'lf(avg)', 'mse',...
         (settings{i}.fs/results{i}.lf*internalsettings.tol_lf(1) + settings{i}.tag.modulation.nfft).^2);
      [errors, errortext] = check_result(i, 2*timings.tsub_s(2), settings{i}.fs/results{i}.lf,...
         errors, errortext, 'lf(std)', 'range', [-Inf, settings{i}.fs/results{i}.lf*internalsettings.tol_lf(2)]);
      [errors, errortext] = check_result(i, timings.t0_s(1)-clock(1)+1, results{i}.t0_s,... % compensate clock shift
         errors, errortext, 'leading t0', 'mse',...
         (results{i}.t0_s*internalsettings.tol_t0 + settings{i}.tag.modulation.nfft/2)^2);
      [errors, errortext] = check_result(i, timings.t0_s(2)+clock(1)-1, results{i}.t0_s,... % compensate clock shift
         errors, errortext, 'trailing t0', 'mse',...
         (results{i}.t0_s*internalsettings.tol_t0 + settings{i}.tag.modulation.nfft/2 +...
         settings{i}.fs/settings{i}.tag.modulation.fclk)^2); % additional tol. for clock shift/jitter
   end
   %     only in case of filtering and we have at least extrap. char. support
   if strcmpi(settings{i}.mode, 'filter') && results{i}.func.pav
      [errors, errortext] = check_result(i, internal.zic, results{i}.zic,...
         errors, errortext, 'Zic', 'relerr', internalsettings.tol_zic);
   end  

   % output (if current loop is last test in block)
   if i == partitioning.indices(index)
      if errors == 0
         disp('         ... passed');
      else
         disp(sprintf('         ... ERRORS: %s', errortext));
         sumoferrors = sumoferrors + 1;  
      end
   end   
end


return
% *******************************************************************************************************
% debugging snippets

% [10*log10(results{i}.pav), 10*log10(internal.pav)]
% [results{i}.pav, internal.pav]
% err_pav = (internal.pav-results{i}.pav)/results{i}.pav
% 
% [10*log10(results{i}.pin), 10*log10(internal.pin)]
% [results{i}.pin, internal.pin]
% err_pin = (internal.pin-results{i}.pin)/results{i}.pin
% 
% [10*log10(results{i}.pic), 10*log10(internal.pic)]
% [results{i}.pic, internal.pic]
% err_pic = (internal.pic-results{i}.pic)/results{i}.pic
% 
% [results{i}.vdda, internal.vdda]
% (internal.vdda-results{i}.vdda)/results{i}.vdda
% 
% disp([(internal.pav-results{i}.pav)/results{i}.pav,...
%    (internal.pin-results{i}.pin)/results{i}.pin,...
%    (internal.pic-results{i}.pic)/results{i}.pic,...
%    (internal.vdda-results{i}.vdda)/results{i}.vdda])
% 
% abs((internal.zic - results{i}.zic) / results{i}.zic)
% 
% [(abs(est_umd(1)) - abs(results{i}.gain(1))) / abs(results{i}.gain(1)), internalsettings.tol_abs(1)]

% clear; close all; clc; pause(0.01);
% version = 'under construction';
% nargin = 1;sumoferrors = 0;errors = 0;errortext = '';
% global globalsettings
% % globalsettings.logging.versions   = 1;
% % globalsettings.logging.exceptions = 1;
% % globalsettings.logging.warnings   = 1;
% % globalsettings.logging.messages   = 1;
% globalsettings.logging.versions   = 0;
% globalsettings.logging.exceptions = 0;
% globalsettings.logging.warnings   = 0;
% globalsettings.logging.messages   = 0;
% 
% % t0
% [settings{i}.tag.modulation.t0*settings{i}.fs, results{i}.t0_s]
% [timings.t0_s(1)-clock(1)+1, timings.t0_s(1)-clock(1)+1-results{i}.t0_s, round(results{i}.t0_s*internalsettings.tol_t0)+settings{i}.tag.modulation.nfft/2]
% % 
% % % magnitude
% [abs(est_umd(1)), abs(results{i}.gain(1)), abs(est_umd(1)) - abs(results{i}.gain(1))]
% [abs(est_mod(1)), abs(results{i}.gain(2)), abs(est_mod(1)) - abs(results{i}.gain(2))]
% % % phase
% % [angle(est_umd(1)), mod(angle(results{i}.gain(1)), pi)]
% % [angle(est_mod(1)), mod(angle(results{i}.gain(2)), pi)]
% % 
% plot
% exp_res1 = nan(size(bitpos));
% exp_res0 = nan(size(bitpos));
% exp_res0(ind_umd) = results{i}.gain(1);
% exp_res1(ind_mod) = results{i}.gain(2);
% figure; hold on;
% plot(magnitude, 'r');
% plot(phase, 'b');
% plot(bitpos, abs(exp_res1), 'ro');
% plot(bitpos, mod(angle(exp_res1), pi), 'bo');
% plot(bitpos, abs(exp_res0), 'rsquare');
% plot(bitpos, mod(angle(exp_res0), pi), 'bsquare');
% legend('abs', 'angle', 'exp:abs (mod)', 'exp:angle (mod)', 'exp:abs (umd)', 'exp:angle (umd)');
% hold off; grid on;
