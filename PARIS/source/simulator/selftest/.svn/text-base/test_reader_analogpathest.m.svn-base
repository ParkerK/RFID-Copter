% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: reader - "analog" path (transmitter, receiver, demodulator) estimator
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
% version = test_reader_analogpathest()
%    Just returns the version number (string).
% sumoferrors = test_reader_analogpathest(sumoferrors)
%    Tests the function reader_analogpathest, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_reader_analogpathest(sumoferrors);
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


function sumoferrors = test_reader_analogpathest(sumoferrors)
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
%     nonrandom tests and symmetry of spectrm
internalsettings.tol_ideal = 1e-15; % MSE(result, expected)
%     random tests
internalsettings.tol_syst = 0.001; % calculated to expected systematic gain (relerr)
internalsettings.tol_rand =  0.05; % calculated to expected random gain (relerr)



% *******************************************************************************************************
% initialization

% output
disp('   = reader_analogpathest **');
disp('     checking returned frequency response at random frequencies (checks are done per bin)');

% load settings and expected results
data = loadmat('results_reader_analogpathest', '     ');
settings     = data.settings;
partitioning = data.partitioning;


% *******************************************************************************************************
% run tests and evaluate

index = 1; % current index in partitioning.indices
for i = 1 : length(settings)
   
   % once for each new test block
   if i == partitioning.indices(index) + 1
      % output
      disp(sprintf('      - %i setup(s) "%s", %i ensemble(s) per setup', partitioning.runs(index),...
         partitioning.names{index}, settings{i}.ens ));
      switch lower(partitioning.names{index})
         case 'ideal'
            disp(sprintf('        checking calculated gains and symmetry MSE < %g', internalsettings.tol_ideal));
         case 'systematic errors'
            disp(sprintf('        checking systematic gain error +/-%g%% and symmetry MSE < %g',...
               internalsettings.tol_syst*100, internalsettings.tol_ideal));
         case 'random errors'
            disp(sprintf('        checking random gain error and symmetry +/-%g%%', internalsettings.tol_rand*100));
         case 'systematic and random errors'
            disp(sprintf('        checking systematic gain error (mean), random gain error (std), and symmetry +/-%g%%',...
               internalsettings.tol_rand*100));
         otherwise
            error('Check partitioning (i=%d).', i);
      end
      % reset errors and errortext
      errors    = 0;
      errortext = '';
      % step index
      index = index + 1;
   end
    
   % ... also have a look at mirrored frequencies to check symmetry
   settings{i}.freq = [-settings{i}.freq; settings{i}.freq];
   
   % calculate expected results
   %     let receiver, demodulator and transmitter report their linear model
   hd_transmitter  = reader_transmitter (settings{i}.reader_transmitter );
   hd_receiver     = reader_receiver    (settings{i}.reader_receiver    );
   hd_demodulation = reader_demodulation(settings{i}.reader_demodulation);
   %     calculate frequency responses
   expres.transmitter  = freqz(hd_transmitter,  settings{i}.fc + settings{i}.freq, settings{i}.fs).'; % at carrier-level
   expres.receiver     = freqz(hd_receiver,     settings{i}.fc + settings{i}.freq, settings{i}.fs).'; % at carrier-level
   expres.demodulation = freqz(hd_demodulation,                  settings{i}.freq, settings{i}.fs).'; % at baseband-level (AAF)
   %     overall frequency response: concatenation = multiplication of spectra
   expres.all = expres.transmitter .* expres.demodulation .* expres.receiver;
   %     get field names
   names_exp = fieldnames(expres);
   
   % call reader_analogpathest and prepare results
   %     for performance purposes: don't call the function settings{i}.ens times with identical settings,
   %     but modify frequency vector instead
   res = reader_analogpathest(repmat(settings{i}.freq, settings{i}.ens, 1), settings{i}.reader_analogpathest,...
      settings{i}.reader_transmitter, settings{i}.reader_receiver, settings{i}.reader_demodulation);
   % check field names of all structs (just to make sure)
   names_res = fieldnames(res);
   if ~all(cellfun(@strcmpi, settings{i}.fieldnames, names_res)) || ~all(cellfun(@strcmpi, settings{i}.fieldnames, names_res))
      error('Fieldnames of returned/expected structs don''t match. Unable to continue the test.');
   end
   % collect ensembles
   for k = 1 : length(settings{i}.fieldnames)
      testres.(settings{i}.fieldnames{k}) = reshape(res.(settings{i}.fieldnames{k}), settings{i}.ens, 2*settings{i}.freq_len);
   end
     
   % checks: all checks done for all fields (and per frequency bin)
   for k = 1 : length(settings{i}.fieldnames)
      % check symmetry (conjugate complex)
      if strcmpi(partitioning.names{index-1}, 'random errors') || strcmpi(partitioning.names{index-1}, 'systematic and random errors')
         %     averaging: allow for larger error
         [errors, errortext] = check_result(i, mean(testres.(settings{i}.fieldnames{k})(:,1:end/2),1),...
            mean(conj(testres.(settings{i}.fieldnames{k})(:,end/2+1:end)),1),...
            errors, errortext, sprintf('symm.%s',settings{i}.fieldnames{k}), 'relerr', internalsettings.tol_rand);
         %     ideal case: should match quite well
      else
         [errors, errortext] = check_result(i, testres.(settings{i}.fieldnames{k})(1:end/2), conj(testres.(settings{i}.fieldnames{k})(end/2+1:end)),...
            errors, errortext, sprintf('symm.%s',settings{i}.fieldnames{k}), 'mse', internalsettings.tol_ideal);
      end
      % no errors
      if strcmpi(partitioning.names{index-1}, 'ideal')
         [errors, errortext] = check_result(i, testres.(settings{i}.fieldnames{k}), expres.(settings{i}.fieldnames{k}),...
            errors, errortext, sprintf('%s',settings{i}.fieldnames{k}), 'mse', internalsettings.tol_ideal);
      end
      % systematic errors (average gain)
      if strcmpi(partitioning.names{index-1}, 'systematic errors')
         [errors, errortext] = check_result(i, testres.(settings{i}.fieldnames{k}),...
            settings{i}.reader_analogpathest.err(1) * expres.(settings{i}.fieldnames{k}),...
            errors, errortext, sprintf('syst.%s',settings{i}.fieldnames{k}), 'relerr', internalsettings.tol_syst);
      elseif strcmpi(partitioning.names{index-1}, 'systematic and random errors')
         [errors, errortext] = check_result(i, mean(testres.(settings{i}.fieldnames{k}), 1),...
            settings{i}.reader_analogpathest.err(1) * expres.(settings{i}.fieldnames{k}),...
            errors, errortext, sprintf('syst.%s',settings{i}.fieldnames{k}), 'relerr', internalsettings.tol_rand);
      end
      % random errors (stddev of gain)
      if strcmpi(partitioning.names{index-1}, 'random errors') || strcmpi(partitioning.names{index-1}, 'systematic and random errors')
         [errors, errortext] = check_result(i,...
            std(testres.(settings{i}.fieldnames{k}) ./ repmat(expres.(settings{i}.fieldnames{k}), settings{i}.ens, 1), 0, 1),...
            settings{i}.reader_analogpathest.err(2) * ones(1, 2*settings{i}.freq_len),...
            errors, errortext, sprintf('rand.%s',settings{i}.fieldnames{k}), 'relerr', internalsettings.tol_rand);
      end         
   end
   
   % clean up
   clear('expres', 'testres');
   
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
%
% % symmetry
% %     rand
% abs((mean(testres.(settings{i}.fieldnames{k})(:,1:end/2),1) -...
%    mean(conj(testres.(settings{i}.fieldnames{k})(:,end/2+1:end)),1)) ./...
%    mean(conj(testres.(settings{i}.fieldnames{k})(:,end/2+1:end))))
% %     ideal
% (mean(testres.(settings{i}.fieldnames{k})(:,1:end/2),1) -...
%    mean(conj(testres.(settings{i}.fieldnames{k})(:,end/2+1:end)),1))
% 
% % ideal
% testres.(settings{i}.fieldnames{k}) - expres.(settings{i}.fieldnames{k})
% 
% % syst. err.
% abs((mean(testres.(settings{i}.fieldnames{k}), 1) -...
%    settings{i}.reader_analogpathest.err(1) * expres.(settings{i}.fieldnames{k}) ) ./...
%    (settings{i}.reader_analogpathest.err(1) * expres.(settings{i}.fieldnames{k})))
% 
% % rand. err
% abs((std(testres.(settings{i}.fieldnames{k}) ./...
%    repmat(expres.(settings{i}.fieldnames{k}), settings{i}.ens, 1), 0, 1) -...
%    settings{i}.reader_analogpathest.err(2)) / settings{i}.reader_analogpathest.err(2))
