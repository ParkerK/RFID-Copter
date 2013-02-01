% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: reader - receiver
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
% version = test_reader_receiver()
%    Just returns the version number (string).
% sumoferrors = test_reader_receiver(sumoferrors)
%    Tests the function reader_receiver, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_reader_receiver(sumoferrors);
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


function sumoferrors = test_reader_receiver(sumoferrors)
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

% specific checks (partitioning has to be checked for changes below!)
% ... check if any of the listed indices matches the current index in partitioning.indices
internalsettings.tests_pwr = [1]; % power splitter tests
internalsettings.tests_bp  = [2]; % bandpass tests

% tolerances
%     general
internalsettings.tol_gain = 0.01; % MSE reported vs. estimated power spectrum (ratio of avg(reported.^2))
%     power splitter tests
internalsettings.tol_ng   = 0.01; % dB noise gain input/output (ratio of expected gain)
%     bandpass filter tests (power spectrum)
internalsettings.tol_pass = 1; % dB passband center (log difference to reported value of power spectrum)


% *******************************************************************************************************
% initialization

% output
disp('   = reader_receiver **');
disp(sprintf('     checking reported vs. estimated power spectrum +/-%g%% of avg(reported.^2)',...
   internalsettings.tol_gain*100));

% load settings and expected results
data = loadmat('results_reader_receiver.mat', '     ');
settings     = data.settings;
results      = data.results;
partitioning = data.partitioning;

% prepare test specific messages
%     check if partitioning names have changed
partitioning.expnames = {'power splitter', 'bandpass + power splitter'};
if length(partitioning.names) ~= length(partitioning.expnames) ||...
      ~all(cellfun(@strcmpi, partitioning.names, partitioning.expnames))
   disp('     ERROR: Partitioning seems to have changed. Skipping tests.'); % see also "check results" part!
end
%     power splitter
partitioning.messages{1} = sprintf('        additional checks: noise gain accuracy +/-%g%%', internalsettings.tol_ng*100);
%     bandpass + power splitter
partitioning.messages{2} = ...
   sprintf('        additional checks: passband center +/-%gdB, stopband gain<=expected', internalsettings.tol_pass);


% *******************************************************************************************************
% run tests and evaluate

index = 1; % current index in partitioning.indices

for i = 1 : length(settings)

   % once for each new test block
   if i == partitioning.indices(index) + 1
      % output
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
      if ~isempty(partitioning.messages{index}); disp(partitioning.messages{index}); end
      % reset errors and errortext
      errors    = 0;
      errortext = '';
      % step index
      index = index + 1;
   end  
   
   % prepare input signal (white noise)
   input = randn(settings{i}.len, 1) * sqrt(settings{i}.var);
   
   % receiver structure
   %     report linear model
   testres{i}.Hd = reader_receiver(settings{i});
   %     work on input signal
   output = reader_receiver(input, settings{i});
   
   % equal vector lengths are a prerequisite for the following calculations => check here
   [errors, errortext] = check_result(i, length(output), settings{i}.len, errors, errortext, 'len', 'equal');
   
   % estimate power spectra of input and I,Q output signals (Welch periodogram estimator)
   %   ... assuming large nfft => unbiased, (except for scaling), influence of window is minimal
   [testres{i}.spec_s, f_s] = pwelch(   input, settings{i}.nfft, settings{i}.nfft/2, settings{i}.nfft, settings{i}.fs);
   [testres{i}.spec_i, f_i] = pwelch(real(output), settings{i}.nfft, settings{i}.nfft/2, settings{i}.nfft, settings{i}.fs);
   [testres{i}.spec_q, f_q] = pwelch(imag(output), settings{i}.nfft, settings{i}.nfft/2, settings{i}.nfft, settings{i}.fs);
   if length(f_s)==length(f_i) && length(f_i)==length(f_q) && all(f_s==f_i) && all(f_i==f_q)
      settings{i}.fvec = f_s;
   end
   %     remove input power spectrum (and most of the estimator bias) => power transfer function
   testres{i}.spec_i = testres{i}.spec_i ./ testres{i}.spec_s;
   testres{i}.spec_q = testres{i}.spec_q ./ testres{i}.spec_s;
   
   % calculate power transfer fcn out of reported linear model for same frequency spacing
   testres{i}.h = abs(freqz(testres{i}.Hd, settings{i}.fvec, settings{i}.fs)).^2;
   
   % compare to expected results
   %     reported vs. estimated spectrum
   %     ... note: comparing the spectra in dB is not a good idea, because log amplifies small differences
   %           => minima of spectrum (cannot be as deep in estimates) will cause large mismatches
   [errors, errortext] = check_result(i, testres{i}.spec_i, testres{i}.h,...
      errors, errortext, 'Hd(I)', 'mse', internalsettings.tol_gain*mean(testres{i}.h.^2));
   [errors, errortext] = check_result(i, testres{i}.spec_q, testres{i}.h,...
      errors, errortext, 'Hd(Q)', 'mse', internalsettings.tol_gain*mean(testres{i}.h.^2));
   %     power splitter tests
   if any(index-1==internalsettings.tests_pwr)
      %     noise gain
      [errors, errortext] = check_result(i, var(real(output)), settings{i}.var*mean(testres{i}.h),... % h is pwr!
         errors, errortext, 'gain(I)', 'relerr', internalsettings.tol_ng*settings{i}.var*mean(testres{i}.h));
      [errors, errortext] = check_result(i, var(imag(output)), settings{i}.var*mean(testres{i}.h),...
         errors, errortext, 'gain(Q)', 'relerr', internalsettings.tol_ng*settings{i}.var*mean(testres{i}.h));
   end
   %     bandpass tests
   if any(index-1==internalsettings.tests_bp)
      %     estimation areas (stopband end/start, log passband center)
      ind.fcut = interp1(settings{i}.fvec, [1:1:length(settings{i}.fvec)], settings{i}.fcut);
      ind.stop = [floor(ind.fcut(1)),  ceil(ind.fcut(2))];
      ind.fc   = interp1(settings{i}.fvec, [1:1:length(settings{i}.fvec)], sqrt(prod(settings{i}.fcut)), 'nearest'); 
      %     stopband attenuation >= set attenuation
      [errors, errortext] = check_result(i,...
         10*log10(max(testres{i}.spec_i([1:ind.stop(1), ind.stop(2):end]))), -settings{i}.att,...
         errors, errortext, 'stopb.att.(I)', 'range', [-Inf, 0]);
      [errors, errortext] = check_result(i,...
         10*log10(max(testres{i}.spec_q([1:ind.stop(1), ind.stop(2):end]))), -settings{i}.att,...
         errors, errortext, 'stopb.att(Q)', 'range', [-Inf, 0]);
      %     passband (log) center attenuation (use reported model: no variance; passb. = 0dB)
      [errors, errortext] = check_result(i, 10*log10(testres{i}.h(ind.fc)), results.att_fc,...
         errors, errortext, 'passb.att', 'range', internalsettings.tol_pass);
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
% figure; hold on;
% plot(settings{i}.fvec/1e6, 10*log10(testres{i}.h), 'b.-')
% plot(settings{i}.fvec/1e6, 10*log10(testres{i}.spec_i), 'r.-');
% plot(settings{i}.fcut/1e6, -ones(size(settings{i}.fcut))*settings{i}.att, 'ko');
% plot(settings{i}.fvec/1e6, 10*log10(abs(testres{i}.h-testres{i}.spec_i)), 'g--');
% plot(settings{i}.fvec/1e6, 10*log10(ones(size(settings{i}.fvec))*mean(abs(testres{i}.h-testres{i}.spec_i))), 'g-');
% hold off; grid on; xlim([min(settings{i}.fvec), max(settings{i}.fvec)]/1e6); ylim([-3*settings{i}.att, 10]);
% legend('expected', 'estimated I', 'fcut', 'error for Hd(I)', 'avg error for Hd(I)');
% ylabel('dB'); xlabel('MHz'); set(gca, 'xScale','log');


