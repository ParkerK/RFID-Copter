% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: reader - demodulation
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
% version = test_reader_demodulation()
%    Just returns the version number (string).
% sumoferrors = test_reader_demodulation(sumoferrors)
%    Tests the function reader_demodulation, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_reader_demodulation(sumoferrors);
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


function sumoferrors = test_reader_demodulation(sumoferrors)
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
%     noise
internalsettings.tol_nvar = 0.15; % ratio of expected noise variance
%     sinusoid and data
%     ...std-dev: anti-aliasing filter as short as possible (performance) => relatively high variance
internalsettings.tol_sym = [0.005, 0.05]; % re/im amplitude (ratio of expected absolute value)
%     ...modulation frequency checked implicitly (symbol bounds calculated out of expected fm)
%     quantization levels
internalsettings.tol_qlev = 0.001; % quantization level tolerance (+/-) to quantization step ratio


% *******************************************************************************************************
% initialization

% store old sumoferrors for later usage
sumoferrors_old = sumoferrors;

% output
disp('   = reader_demodulation **');

% load settings and expected results
data = loadmat('results_reader_demodulation.mat', '     ');
settings     = data.settings;
oscsettings  = data.oscsettings;
partitioning = data.partitioning;


% *******************************************************************************************************
% run tests and evaluate

index = 1; % current index in partitioning.indices
runtime_warn.long_impres = []; % indices of test runs with demod filter impres > symbol duration

for i = 1 : length(settings)

   % once for each new test block
   if i == partitioning.indices(index) + 1
      % output
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
      switch lower(partitioning.names{index})
         case 'noise'
            disp(sprintf('        checking noise gain +/-%g%%', internalsettings.tol_nvar*100));
         case{'sinusoid'};
            disp(sprintf('        checking I/Q signals avg +/-%g%%, std<=%g%%',...
            internalsettings.tol_sym(1)*100, internalsettings.tol_sym(2)*100));
         case{'data'};
            disp(sprintf('        checking I/Q signals avg +/-%g%%, std<=%g%% for each symbol',...
            internalsettings.tol_sym(1)*100, internalsettings.tol_sym(2)*100));
         case 'quantization'
            disp(sprintf('        checking quantization steps +/-%g%% of amplitude resolution', internalsettings.tol_qlev*100));
      end
      % reset errors and errortext
      errors    = 0;
      errortext = '';
      % step index
      index = index + 1;
   end
   
   % prepare input vector (if-structures only to speed up the process)
   input = zeros(settings{i}.veclen_s, 1);
   %     carrier?
   if settings{i}.mode.carrier
      input = input + sqrt(2*settings{i}.varc)*exp(complex(0,2*pi*settings{i}.f0 * [0:1:settings{i}.veclen_s-1]'/settings{i}.fs));
   end
   %     data? (complex modulation)
   if settings{i}.mode.data
      % modulation signal (@ fm: settings{i}.data)
      modsignal = zeros(settings{i}.veclen_s, 1);
      for symb = 1 : settings{i}.symb
         ind = 1+round((symb-1)*settings{i}.fs/settings{i}.fm:symb*settings{i}.fs/settings{i}.fm);
         modsignal(ind) = settings{i}.data(symb) * ones(length(ind), 1); % complex mod
      end
      % modulate
      input(1:length(modsignal)) = input(1:length(modsignal)) .* modsignal;
   end
   %     noise?
   if settings{i}.mode.noise
      input = input + randn(settings{i}.veclen_s, 1) * sqrt(settings{i}.varn);
   end
         
   % demodulate
   [output, internal] = reader_demodulation(complex(real(input), real(input)), settings{i}, oscsettings{i});
      
   % get linear baseband model of reader_demodulation (= anti-aliasing filter)
   Hd = reader_demodulation(settings{i});
   %     calculate gain factor for f0=dc and f0+fm=fm (attention: freqz needs a frequency VECTOR)
   testres{i}.h = freqz(Hd, [0, settings{i}.fm], settings{i}.fs);
   %     impulse response (let impres choose the needed length) plus sampled version
   testres{i}.impres    = impz(Hd);
   testres{i}.impres_rs = testres{i}.impres(round(1:settings{i}.fs/settings{i}.frs:length(testres{i}.impres)-1));
   %     noise gain
   testres{i}.ng = sum(testres{i}.impres.^2);
   
   % check for impres > symbol duration and correct if necessary
   if length(testres{i}.impres_rs) > round(settings{i}.frs/settings{i}.fm)
      runtime_warn.long_impres = [runtime_warn.long_impres, i];
   end
     
   % expected modsignal (sampled, but unquantized)
   if ~settings{i}.mode.data
      results{i}.modsignal_rs = ones(size(output)) * sqrt(settings{i}.varc) * testres{i}.h(1); % carrier (mapped to DC)
   else
      results{i}.modsignal_rs = modsignal(round(1:settings{i}.fs/settings{i}.frs:length(modsignal)-1)) * sqrt(settings{i}.varc);
   end    
     
   % (re-)initialize range arrays
   for j = min(settings{i}.ind_data) : max(settings{i}.ind_data) % for all possible constellations
      symbrange{j} = [];
   end
    
   % select ranges for estimators
   %     data ... symbol periods do not have to be equal (e.g. frs/fm not integer) => not vectorized
   if settings{i}.mode.data
      % calculate range for each symbol
      for j = 1 : length(settings{i}.data) % end of filter impres to end of symbol
         % expected start/end of steady symbol (ignored transient phase)
         symbbord = round([j-1, j]*settings{i}.frs/settings{i}.fm + [length(testres{i}.impres_rs), -1]);
         % assure operability in case of too long impulse responses (check for warning above)
         if symbbord(1) > symbbord(2)
            symbbord(1) = symbbord(1);
         end
         % add range to corresponding symbol
         symbrange{settings{i}.ind_data(j)} = [symbrange{settings{i}.ind_data(j)}, symbbord(1):symbbord(2)];
      end
   %     no data => use everything from end of transient phase (impres) to start of zero-padding
   %     ... in this case, symbrange{1} does not correspond to a mod constellation!
   %     ... make sure this works for truncated and zero-padded signals as well
   else
      symbrange{1} = 1+length(testres{i}.impres_rs):floor(...
         min(settings{i}.veclen_s*settings{i}.frs/settings{i}.fs, length(output)));
   end
     
   % estimate
   %     noise
   if strcmpi(partitioning.names{index-1}, 'noise')
      testres{i}.varn = [var(real(output(symbrange{1}))); var(imag(output(symbrange{1})))];
   else
      testres{i}.varn = NaN;
   end
   %     carrier only: average I/Q channels (only one "symbol"=1)
   %     data: average I/Q channels for all symbols (vectors might be empty)
   if strcmpi(partitioning.names{index-1}, 'sinusoid') || strcmpi(partitioning.names{index-1}, 'data')
      for j = 1 : length(symbrange)
         est{j} = [complex(mean(real(output(symbrange{j}))), mean(imag(output(symbrange{j})))),...
                   complex( std(real(output(symbrange{j}))),  std(imag(output(symbrange{j}))))];
      end
   end
   %     quantization: create histogram; position of quantization steps
   if strcmpi(partitioning.names{index-1}, 'quantization')      
      % create normalized class edges for histogram
      % ... will result in the following classes (ql: quantization level, Nq=2^settings{i}.q): 
      %     OUT:<ql0 | IN:ql0 | OUT:ql0-ql1 | IN:ql1 | ... | IN:qlNq | OUT:>qlNq
      % ... thus: all odd classes (histc) must be zero (OUT)
      testres{i}.q_classes.norm = [0 : 1 : 2^settings{i}.q-1] - 2^(settings{i}.q-1) + 0.5;
      testres{i}.q_edges.norm = ...
         [-Inf, sort([testres{i}.q_classes.norm - internalsettings.tol_qlev*2^-settings{i}.q,...
         testres{i}.q_classes.norm + internalsettings.tol_qlev*2^-settings{i}.q]), Inf];
      % scale edges to bounds of I (re) and Q (im) channel of demodulated signal
      testres{i}.q_edges.re = testres{i}.q_edges.norm / internal.agc_re;
      testres{i}.q_edges.im = testres{i}.q_edges.norm / internal.agc_im;
      % create histograms
      testres{i}.q_hist.re = histc(real(output), testres{i}.q_edges.re);
      testres{i}.q_hist.im = histc(imag(output), testres{i}.q_edges.im);
      %      combine the last two entries to get <= Inf and delete last entry
      testres{i}.q_hist.re(end-1) = testres{i}.q_hist.re(end-1) +  testres{i}.q_hist.re(end);
      testres{i}.q_hist.im(end-1) = testres{i}.q_hist.im(end-1) +  testres{i}.q_hist.im(end);
      testres{i}.q_hist.re(end) = [];
      testres{i}.q_hist.im(end) = [];
      % check histograms
      %      "forbidden" zone between levels (OUT)
      testres{i}.q_hist.out_re = sum(testres{i}.q_hist.re(1:2:end));
      testres{i}.q_hist.out_im = sum(testres{i}.q_hist.im(1:2:end));
      %      within quantization levels (IN)
      testres{i}.q_hist.in_re = sum(testres{i}.q_hist.re(2:2:end-1));
      testres{i}.q_hist.in_im = sum(testres{i}.q_hist.im(2:2:end-1));
   end
            
   % compare to expected results
   %     vector length (if not automatically set); equality I/Q channel signal lengths above
   if ~isnan(settings{i}.len_rs)
      [errors, errortext] = check_result(i, length(output), settings{i}.len_rs,...
         errors, errortext, 'len', 'equal');
   end
   %     noise: noise gain
   if strcmpi(partitioning.names{index-1}, 'noise')
      [errors, errortext] = check_result(i, testres{i}.varn(1)/testres{i}.ng, settings{i}.varn,...
         errors, errortext, 'varn,re', 'relerr', internalsettings.tol_nvar);
      [errors, errortext] = check_result(i, testres{i}.varn(2)/testres{i}.ng, settings{i}.varn,...
         errors, errortext, 'varn,im', 'relerr', internalsettings.tol_nvar);
   end
   %     carrier only (would be identical to symbol=1)
   if strcmpi(partitioning.names{index-1}, 'sinusoid')
      %     averages (real and imaginary part; ATTENTION: this is not equal to relerr)
      [errors, errortext] = check_result(i, real(est{1}(1)), mean(real(results{i}.modsignal_rs(symbrange{1}))),...
         errors, errortext, 'avg,re', 'mse',...
         (internalsettings.tol_sym(1)*mean(abs(results{i}.modsignal_rs(symbrange{1})))).^2);
      [errors, errortext] = check_result(i, imag(est{1}(1)), mean(imag(results{i}.modsignal_rs(symbrange{1}))),...
         errors, errortext, 'avg,im', 'mse',...
         (internalsettings.tol_sym(1)*mean(abs(results{i}.modsignal_rs(symbrange{1})))).^2);
      %     standard deviations (real and imaginary part)
      [errors, errortext] = check_result(i, real(est{1}(2)), 0,...
         errors, errortext, 'std,re', 'range',...
         internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{1}))));
      [errors, errortext] = check_result(i, imag(est{1}(2)), 0,...
         errors, errortext, 'std,im', 'range',...
         internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{1}))));
   end
   %     data: symbol j
   if strcmpi(partitioning.names{index-1}, 'data')
      for j = 1 : length(symbrange)
         %     averages (real and imaginary part; ATTENTION: this is not equal to relerr)
         [errors, errortext] = check_result(i, real(est{j}(1)), mean(real(results{i}.modsignal_rs(symbrange{j}))),...
            errors, errortext, sprintf('avg,re(sym%i)',j), 'mse',...
            (internalsettings.tol_sym(1)*mean(abs(results{i}.modsignal_rs(symbrange{j})))).^2);
         [errors, errortext] = check_result(i, imag(est{j}(1)), mean(imag(results{i}.modsignal_rs(symbrange{j}))),...
            errors, errortext, sprintf('avg,im(sym%i)',j), 'mse',...
            (internalsettings.tol_sym(1)*mean(abs(results{i}.modsignal_rs(symbrange{j})))).^2);
         %     standard deviations (real and imaginary part)
         [errors, errortext] = check_result(i, real(est{j}(2)), 0,...
            errors, errortext, sprintf('std,re(sym%i)',j), 'range',...
            internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{j}))));
         [errors, errortext] = check_result(i, imag(est{j}(2)), 0,...
            errors, errortext, sprintf('std,im(sym%i)',j), 'range',...
            internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{j}))));
      end
   end
   %     quantization: quantization levels
   if strcmpi(partitioning.names{index-1}, 'quantization')
      %     no values in between quantization steps allowed
      [errors, errortext] = check_result(i, testres{i}.q_hist.out_re, 0,...
         errors, errortext, 'val,re', 'equal');
      [errors, errortext] = check_result(i, testres{i}.q_hist.out_re, 0,...
         errors, errortext, 'val,im', 'equal'); 
      %     double-check histogram results
      [errors, errortext] = check_result(i, testres{i}.q_hist.in_re+testres{i}.q_hist.out_re, length(output),...
         errors, errortext, 'TEST ERROR: HIST,re', 'equal');
      [errors, errortext] = check_result(i, testres{i}.q_hist.in_im+testres{i}.q_hist.out_im, length(output),...
         errors, errortext, 'TEST ERROR: HIST,im', 'equal');
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

% Additional message: impres > symbol duration warning
if ~isempty(runtime_warn.long_impres) && sumoferrors > sumoferrors_old
   disp('   Problems might have been caused by too long demod filter impulse response for indices:');
   indices = sprintf('%i,', runtime_warn.long_impres);
   disp(sprintf('      %s', indices(1:end-1)));
end


return
% *******************************************************************************************************
% debugging snippets

% % clear; close all; clc; version = 'UC'; nargin = 1; sumoferrors = 0;
% 
% close all;
% 
% % I-channel plot (does NOT work for noise only!)
% figure; subplot(2,1,1); hold on;
% plot(real(results{i}.modsignal_rs), 'b'); % expected
% plot(real(output), 'r'); % calculated result
% plot(testres{i}.impres_rs/max(testres{i}.impres_rs)*max(real(results{i}.modsignal_rs)), 'g'); % Hd impulse response
% for j = 1 : length(symbrange); plot(symbrange{j}, real(est{j}(1))*ones(size(symbrange{j})), 'k.'); end % estimation areas
% hold off; grid on; axis tight;
% title('I-CHANNEL'); xlabel('samples @ frs'); legend('expected', 'calculated', 'Impres', 'est. areas');
%    subplot(2,1,2); hold on;
% plot(real(output)-real(results{i}.modsignal_rs), 'b');
% plot( ones(size(output)) * internalsettings.tol_sym(1)*mean(abs(results{i}.modsignal_rs(symbrange{j}))), 'g-'); % bound avg
% plot( ones(size(output)) * internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{j}))), 'g--'); % bound std
% for j = 1 : length(symbrange); plot(symbrange{j}, zeros(size(symbrange{j})), 'k.'); end % estimation areas
% plot(-ones(size(output)) * internalsettings.tol_sym(1)*mean(abs(results{i}.modsignal_rs(symbrange{j}))), 'g-'); % bound avg
% plot(-ones(size(output)) * internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{j}))), 'g--'); % bound std
% hold off; grid on; axis tight;
% title('ERROR I-CHANNEL'); xlabel('samples @ frs'); legend('error', 'bound avg', 'bound std', 'est. areas');
% ylim([-1.2,1.2]*internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{j}))));
% 
% % Q-channel plot (does NOT work for noise only!)
% figure; subplot(2,1,1); hold on;
% plot(imag(results{i}.modsignal_rs), 'b'); % expected
% plot(imag(output), 'r'); % calculated result
% plot(testres{i}.impres_rs/max(testres{i}.impres_rs)*max(imag(results{i}.modsignal_rs)), 'g'); % Hd impulse response
% for j = 1 : length(symbrange); plot(symbrange{j}, imag(est{j}(1))*ones(size(symbrange{j})), 'k.'); end % estimation areas
% hold off; grid on; axis tight;
% title('Q-CHANNEL'); xlabel('samples @ frs'); legend('expected', 'calculated', 'Impres', 'est. areas');
%    subplot(2,1,2); hold on;
% plot(imag(output)-imag(results{i}.modsignal_rs), 'b'); % error
% plot( ones(size(output)) * internalsettings.tol_sym(1)*mean(abs(results{i}.modsignal_rs(symbrange{j}))), 'g-'); % bound avg
% plot( ones(size(output)) * internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{j}))), 'g--'); % bound std
% for j = 1 : length(symbrange); plot(symbrange{j}, zeros(size(symbrange{j})), 'k.'); end % estimation areas
% plot(-ones(size(output)) * internalsettings.tol_sym(1)*mean(abs(results{i}.modsignal_rs(symbrange{j}))), 'g-'); % bound avg
% plot(-ones(size(output)) * internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{j}))), 'g--'); % bound std
% hold off; grid on; axis tight; legend('expected');
% title('ERROR Q-CHANNEL'); xlabel('samples @ frs'); legend('error', 'bound avg', 'bound std', 'est. areas');
% ylim([-1.2,1.2]*internalsettings.tol_sym(2)*mean(abs(results{i}.modsignal_rs(symbrange{j}))));
% 
% % find intersections with 50% threshold for timing detection
% am50 = [min(real(results{i}.modsignal_rs))+max(real(results{i}.modsignal_rs));...
%    min(imag(results{i}.modsignal_rs))+max(imag(results{i}.modsignal_rs))] / 2; % [I; Q]
% testres{i}.int50i_rs = findzeros([-1; real(output)-am50(1); -1])'; % | not necessarily
% testres{i}.int50q_rs = findzeros([-1; imag(output)-am50(1); -1])'; % | of equal length!

%       figure; hold on;
%       plot(real(output), 'b.-');
%       plot([1, length(output)], [testres{i}.q_edges.re', testres{i}.q_edges.re']);
%       hold off; grid on;
%
%       figure; hold on;
%       plot(imag(output), 'b.-');
%       plot([1, length(output)], [testres{i}.q_edges.im', testres{i}.q_edges.im']);
%       hold off; grid on;
