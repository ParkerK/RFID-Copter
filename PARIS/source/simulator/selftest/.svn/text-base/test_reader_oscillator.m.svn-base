% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: reader - oscillator
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
% version = test_reader_oscillator()
%    Just returns the version number (string).
% sumoferrors = test_reader_oscillator(sumoferrors)
%    Tests the function reader_oscillator, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_reader_oscillator(sumoferrors);
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
% - improve performance of cleancarrier and phasedev subfunctions
% ? time-variant checks
%
% *******************************************************************************************************


function sumoferrors = test_reader_oscillator(sumoferrors)
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
%     clean carrier
internalsettings.tol_a = 1e-3; % relative amplitude error
internalsettings.tol_f = 1e-4; % relative frequency error
internalsettings.tol_p = 1e-2; % degree absolute phase error
%     instabilities, noise
internalsettings.tol_astddev = 2e-2; % relative error, amplitude instability
internalsettings.tol_fstddev = 2e-2; % relative error, frequency instability
internalsettings.tol_snr     = 5e-1; % dB abolute error, carrier to noise ratio


% *******************************************************************************************************
% initialize

% output
disp('   = reader_oscillator **');

% load settings and expected results
data = loadmat('results_reader_oscillator.mat', '     ');
settings     = data.settings;
results      = data.results;
partitioning = data.partitioning;


% *******************************************************************************************************
% run tests and evaluate

for i = 1 : length(settings)
   
   % partitioning and output
   %     find index for partition struct (partitioning.indices is "end-of-block")
   for j = 1 : length(partitioning.indices)
      if i <= partitioning.indices(j)
         index = j;
         break;
      end
   end
   
   % display test mode / reset errors for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
      switch lower(partitioning.names{index})
         case 'clean carrier'
            disp(sprintf('        checking frequency +/-%g%%, phase +/-%gdeg, and amplitude +/-%g%%',...
               internalsettings.tol_f*100, internalsettings.tol_p, internalsettings.tol_a*100));
         case 'amplitude instability'
            disp(sprintf('        checking amplitude instability +/-%g%%',...
               internalsettings.tol_astddev*100));
         case 'frequency instability'
            disp(sprintf('        checking frequency instability +/-%g%%',...
               internalsettings.tol_fstddev*100));
         case 'additive noise'
            disp(sprintf('        checking SNR +/-%gdB', internalsettings.tol_snr));
         otherwise
            error('Unsupported test "%s".', lower(partitioning.names{index}));
      end
      % reset errors
      errors    = 0;
      errortext = '';
   end
   
   % call reader_oscillator
   carrier = reader_oscillator(settings{i});
   
   % checks
   switch lower(partitioning.names{index})
      case 'clean carrier'
         % get (average) frequency out of spectrum
         nfft = 2^nextpow2(settings{i}.fs * 2 / (internalsettings.tol_f * settings{i}.fcenter));
         X = fft(carrier, nfft);
         [dummy, ind_max] = max(abs(X));
         testresults.f = (ind_max-1) / nfft * settings{i}.fs;
         % get amplitude and phase (est_sinusoid is relative to a cosine => modify phase)
         testresults.ap_re = est_sinusoid(real(carrier),...
            struct('n',length(carrier), 'ol', 0, 'f0',settings{i}.fcenter/settings{i}.fs));
         testresults.ap_re.p = testresults.ap_re.p + 90;
         if strcmpi(settings{i}.mode, 'exp')
            testresults.ap_im = est_sinusoid(imag(carrier),...
               struct('n',length(carrier), 'ol', 0, 'f0',settings{i}.fcenter/settings{i}.fs));
            testresults.ap_im.p = testresults.ap_im.p + 90;
         end
         % checks
         [errors, errortext] = check_result(i, testresults.f, settings{i}.fcenter,...
            errors, errortext, 'freq', 'relerr', internalsettings.tol_f);
         [errors, errortext] = check_result(i, testresults.ap_re.a, results{i}.ampl,...
            errors, errortext, 'ampl(re)', 'relerr', internalsettings.tol_a);
         [errors, errortext] = check_result(i, testresults.ap_re.p, results{i}.phase_re,...
            errors, errortext, 'phase(re)', 'range', [-1,1]*internalsettings.tol_p);
         if strcmpi(settings{i}.mode, 'exp')
            [errors, errortext] = check_result(i, testresults.ap_re.a, results{i}.ampl,...
               errors, errortext, 'ampl(im)', 'relerr', internalsettings.tol_a);
            [errors, errortext] = check_result(i, testresults.ap_re.p, results{i}.phase_re,...
               errors, errortext, 'phase(im)', 'range', [-1,1]*internalsettings.tol_p);
         end
   
      case 'amplitude instability'
         ccarrier = cleancarrier(length(carrier), settings{i}, results{i});
         testresults.astddev = nanstd((carrier - ccarrier) ./ ccarrier);
         [errors, errortext] = check_result(i, testresults.astddev, settings{i}.astddev,...
            errors, errortext, 'astddev', 'relerr', internalsettings.tol_astddev);

      case 'frequency instability'
         ccarrier = cleancarrier(length(carrier), settings{i}, results{i});
         testresults.pdev    = phasedev(carrier, ccarrier, settings{i});
         testresults.fstddev = std(testresults.pdev)  * settings{i}.fcenter/(2*pi); % Hz
         [errors, errortext] = check_result(i, testresults.fstddev, settings{i}.fstddev,...
            errors, errortext, 'fstddev', 'relerr', internalsettings.tol_fstddev);
         
      case 'additive noise'
         testresults.snr = 10*log10( var(carrier) / ...
            var(carrier - cleancarrier(length(carrier), settings{i}, results{i})) );
         [errors, errortext] = check_result(i, testresults.snr, settings{i}.snr,...
            errors, errortext, 'SNR', 'range', [-1,1]*internalsettings.tol_snr);
         
      otherwise
         error('Unsupported test "%s".', lower(partitioning.names{index}));
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

if sumoferrors > 0
   return
end

end



% *******************************************************************************************************
% *******************************************************************************************************
% create a clean carrier (for amplitude/phase noise tests)
%    len_s   length of signal in samples @ set.fs
%    set     current settings struct
%    res     current results struct
function signal = cleancarrier(len_s, set, res)
switch lower(set.mode)
   case 'cos'
      signal = res.ampl * cos(2*pi * set.fcenter/set.fs * [0:len_s-1]');
   case 'sin'
      signal = res.ampl * sin(2*pi * set.fcenter/set.fs * [0:len_s-1]');
   case 'exp'
      signal = res.ampl * exp(complex(0, 2*pi * set.fcenter/set.fs * [0:len_s-1]'));
   otherwise
      error('Unsupported mode="%s".', lower(set.mode));
end
end


% *******************************************************************************************************
% *******************************************************************************************************
% deviation of phase to clean carrier (for phase noise tests)
%    carrier        signal created by reader_oscillator
%    cleancarrier   signal without phase noise
%    set            current settings struct
function pdev = phasedev(carrier, cleancarrier, set)
switch lower(set.mode)
   case 'cos'
      pdev = acos(carrier/sqrt(2)) - acos(cleancarrier/sqrt(2));
   case 'sin'
      pdev = asin(carrier/sqrt(2)) - asin(cleancarrier/sqrt(2));
   case 'exp'
      pdev = imag(log(carrier./cleancarrier));
   otherwise
      error('Unsupported mode="%s".', lower(set.mode));
end
end

