% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: reader - transmitter
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
% version = test_reader_transmitter()
%    Just returns the version number (string).
% sumoferrors = test_reader_transmitter(sumoferrors)
%    Tests the function reader_transmitter, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_reader_transmitter(sumoferrors);
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
%
% *******************************************************************************************************


function sumoferrors = test_reader_transmitter(sumoferrors)
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

% tolerances
%     general
internalsettings.tol_gain = 0.05; % MSE reported vs. estimated power spectrum (ratio of avg(reported.^2))
%     power splitter tests
internalsettings.tol_ng   = 0.01; % dB noise gain input/output (ratio of expected gain)
%     bandpass tests
internalsettings.tol_att  = 0.01; % dB noise gain input/output (ratio of expected gain)


% *******************************************************************************************************
% initialization

% output
disp('   = reader_transmitter **');
disp(sprintf('     checking reported vs. estimated power spectrum +/-%g%% of avg(reported.^2)',...
   internalsettings.tol_gain*100));

% load settings and expected results
data = loadmat('results_reader_transmitter', '     ');
settings     = data.settings;
partitioning = data.partitioning;


% *******************************************************************************************************
% run tests and evaluate

index = 1; % current index in partitioning.indices
for i = 1 : length(settings)

   % once for each new test block
   if i == partitioning.indices(index) + 1
      % output
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
      switch lower(partitioning.names{index})
         case {'linear amp', 'pwramp'}
            disp(sprintf('        additional check: noise gain accuracy +/-%g%%', internalsettings.tol_ng*100));
         case {'linear amp + bandpass', 'pwramp + bandpass'}
            disp(sprintf('        additional check: stopband gain<=expected (%g dB tolerance)', internalsettings.tol_att));
         otherwise
            error('Check partitioning (i=%d).', i);
      end
      if strcmpi(partitioning.names{index}, 'pwramp') || strcmpi(partitioning.names{index}, 'pwramp + bandpass') 
         disp('        Using pre-warping by expected pwramp characteristic to deal with nonlinearity.');
      end
      % reset errors and errortext
      errors    = 0;
      errortext = '';
      % step index
      index = index + 1;
   end  
   
   % prepare input signal (white noise)
   input = randn(settings{i}.len, 1) * sqrt(settings{i}.var);
   %     pre-warp to counter effects of nonlinear power amplifier if enabled (pwramp is before bandpass)
   if settings{i}.en_nl
      % bounds of signal
      bounds = [min(input), max(input)];
      % load characteristic
      pwramp = loadmat(settings{i}.nl_char);
      %     check for symmetry
      if abs(mean(pwramp.nl_x)) > eps || abs(mean(pwramp.nl_y)) > eps
         criterr('Characteristic of power-amplifier has to be symmetrical around zero.');
      end
      % invert characteristic nl_y(nl_x) and scale such that inv_x=-1...1
      scale = 1 / (max(pwramp.nl_y)-min(pwramp.nl_y)) * (max(pwramp.nl_x)-min(pwramp.nl_x));
      pwramp.inv_x = pwramp.nl_y * scale;
      pwramp.inv_y = pwramp.nl_x * scale;
      %     warp limited input swing
      pwramp.inv_lim = 2*interp1(pwramp.nl_x, pwramp.nl_y, settings{i}.nl_lim, 'linear', 'extrap') / ...
         (max(pwramp.nl_y)-min(pwramp.nl_y));
      % scale characteristic to input signal (symmetrical)
      %     for x: -1 ... 1 => 0 ... 1 plus apply scaling by inv_lim (gain) and nl_lim
      pwramp.inv_x = (pwramp.inv_x/pwramp.inv_lim                   +1)/2;
      pwramp.inv_y = (pwramp.inv_y/pwramp.inv_lim/settings{i}.nl_lim+1)/2;
      %     for x: "0 ... 1 => min(signal) ... max(signal) "
      pwramp.inv_x = pwramp.inv_x * diff(bounds) + bounds(1);
      pwramp.inv_y = pwramp.inv_y * diff(bounds) + bounds(1);
      % pre-warp (inverse interp1 of reader_transmitter)
      input_tx = interp1(pwramp.inv_x, pwramp.inv_y, input, 'linear', 'extrap');
   else
      input_tx = input;
   end
        
   % transmitter structure
   %     report linear model
   testres{i}.Hd = reader_transmitter(settings{i});
   %     work on input signal
   output = reader_transmitter(input_tx, settings{i});
      
   % equal vector lengths are a prerequisite for the following calculations => check here
   [errors, errortext] = check_result(i, length(output), settings{i}.len, errors, errortext, 'len', 'equal');
   
   % estimate power spectra of input and I,Q output signals (Welch periodogram estimator)
   %   ... assuming large nfft => unbiased, (except for scaling), influence of window is minimal
   [testres{i}.spec_i, f_i] = pwelch( input, settings{i}.nfft, settings{i}.nfft/2, settings{i}.nfft, settings{i}.fs);
   [testres{i}.spec_o, f_o] = pwelch(output, settings{i}.nfft, settings{i}.nfft/2, settings{i}.nfft, settings{i}.fs);
   if length(f_i)==length(f_o) && all(f_i==f_o)
      settings{i}.fvec = f_i;
   else
      error('Frequency vectors created by pwelch are not equal (i=%d).', i);
   end
   %     remove input power spectrum (and most of the estimator bias) => power transfer function
   testres{i}.spec_o = testres{i}.spec_o ./ testres{i}.spec_i;
   
   % calculate power transfer fcn out of reported linear model for same frequency spacing
   % ... variance of input is removed by final power normalization, but not known to reader_transmitter
   %     beforehand => remove here
   % ... nonlinearity is cancelled before the filter => the system should be linear
   testres{i}.h_pwr = abs(freqz(testres{i}.Hd, settings{i}.fvec, settings{i}.fs)).^2 / settings{i}.var;
   
   % compare to expected results
   %     reported vs. estimated spectrum
   %     ... note: comparing the spectra in dB is not a good idea, because log amplifies small differences
   %           => minima of spectrum (cannot be as deep in estimates) will cause large mismatches
   [errors, errortext] = check_result(i, testres{i}.spec_o, testres{i}.h_pwr,...
      errors, errortext, 'Hd', 'mse', internalsettings.tol_gain*mean(testres{i}.h_pwr.^2));
   switch lower(partitioning.names{index-1})
      case {'linear amp', 'pwramp'}
         %     noise gain
         [errors, errortext] = check_result(i, var(output), settings{i}.var*mean(testres{i}.h_pwr),... % h is pwr!
            errors, errortext, 'gain', 'relerr', internalsettings.tol_ng*settings{i}.var*mean(testres{i}.h_pwr));
      case {'linear amp + bandpass', 'pwramp + bandpass'}
         %     estimation areas (stopband end/start, log passband center)
         ind.fcut = interp1(settings{i}.fvec, [1:1:length(settings{i}.fvec)], settings{i}.bp_fcut);
         ind.stop = [floor(ind.fcut(1)),  ceil(ind.fcut(2))];
         %     stopband attenuation >= set attenuation (use reported model: no variance)
         %     (still there is some normalization in between => add small tolerance)
         [errors, errortext] = check_result(i,...
            10*log10(max(testres{i}.h_pwr((ind.stop(1)+1:ind.stop(2)-1)))) - ...
            10*log10(max(testres{i}.h_pwr([1:ind.stop(1), ind.stop(2):end]))), settings{i}.bp_att,...
            errors, errortext, 'stopb.att.', 'range', [-internalsettings.tol_att, Inf]);
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

% clear; close all; clc; version = 'UC'; nargin = 1; sumoferrors = 0;

% close all;
% 
% % averaged spectra
% bp_bordy = [-1.05*settings{i}.bp_att, 1.05*max(10*log10(testres{i}.spec_o))];
% figure; hold on;
% %     expected
% plot(settings{i}.fvec/1e6, 10*log10(testres{i}.h_pwr), 'b.-')
% %     estimated
% plot(settings{i}.fvec/1e6, 10*log10(testres{i}.spec_o), 'r.-');
% %     fcut
% plot(ones(2,1)*settings{i}.bp_fcut(1)/1e6, bp_bordy, 'k-');
% %     stopband gain
% plot(settings{i}.fvec/1e6,...
%    ones(size(settings{i}.fvec))*10*log10(max(testres{i}.h_pwr((ind.stop(1)+1:ind.stop(2)-1)))), 'k-');
% %     error Hd
% plot(settings{i}.fvec/1e6, 10*log10(abs(testres{i}.h_pwr-testres{i}.spec_o)), 'g--');
% %     avg error Hd
% plot(settings{i}.fvec/1e6,...
%    10*log10(ones(size(settings{i}.fvec))*mean(abs(testres{i}.h_pwr-testres{i}.spec_o))), 'g-');
% %     fcut and stopband gain II
% plot(ones(2,1)*settings{i}.bp_fcut(2)/1e6, bp_bordy, 'k-');
% plot(settings{i}.fvec/1e6,...
%    ones(size(settings{i}.fvec))*10*log10(max(testres{i}.h_pwr((ind.stop(1)+1:ind.stop(2)-1))))-...
%    settings{i}.bp_att, 'k-');
% %     ...
% hold off; grid on; xlim([min(settings{i}.fvec), max(settings{i}.fvec)]/1e6);
% ylim([-1.1*settings{i}.bp_att, 1.1*max(10*log10(testres{i}.spec_o))]);
% legend('expected', 'estimated', 'fcut', 'stopb.gain', 'error for Hd', 'avg error for Hd');
% ylabel('dB'); xlabel('MHz'); set(gca, 'xScale','log');
% 
%   
% [mean((testres{i}.spec_o-testres{i}.h_pwr).^2),internalsettings.tol_gain*mean(testres{i}.h_pwr.^2)]
% mean((testres{i}.spec_o-testres{i}.h_pwr).^2)/(internalsettings.tol_gain*mean(testres{i}.h_pwr.^2))

