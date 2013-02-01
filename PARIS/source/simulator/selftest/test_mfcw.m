% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% (quick) selftest: ranging - multi-frequency continuous-wave (MFCW) ranging
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
% version = test_mfcw()
%    Just returns the version number (string).
% sumoferrors = test_mfcw(sumoferrors)
%    Tests the functions mfcw_addseccarrier, mfcw_compsel, and mfcw_calcdist, displays the results in  
%    the command window and returns the overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_mfcw(sumoferrors);
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
% ? complete tests such that all parameters are checked (gammai, c_ord, ...)
%
% *******************************************************************************************************


function sumoferrors = test_mfcw(sumoferrors)
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

% anti-aliasing filter (filtfilt)
internalsettings.aaf_ord =  4; % filter order
internalsettings.aaf_att = 80; % dB attenuation

% range for distance estimates (skip filter transients)
internalsettings.est_ind = [1/3, 2/3]; % [min, max] 0...1

% result tolerances
internalsettings.tol = 0.05; % tolerance for distance (relative error)


% *******************************************************************************************************
% initialization

% temporarily switch off MATLAB:nearlySingularMatrix (some setting might result in bad matrix scaling)
oldwarn = warning('off', 'MATLAB:nearlySingularMatrix');

% output
disp('   = multi-frequency-continuous-wave ranging (mfcw_addseccarrier, mfcw_compsel, and mfcw_calcdist) ***');
disp(sprintf('     checking estimated distance +/-%g%%', internalsettings.tol*100));

% load settings and expected results
data = loadmat('results_mfcw', '     ');
settings     = data.settings;
partitioning = data.partitioning;


% *******************************************************************************************************
% run tests and evaluate

for i = 1 : length(partitioning.names)
   
   % once for each new test block
   %     output
   disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(i), partitioning.names{i}));
   %     reset errors and errortext
   errors    = 0;
   errortext = '';
   
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)
      
      % signal path
      %     create time vector and base carriers
      t = [0 : 1 / settings{j}.fs : (settings{j}.len_s-1) / settings{j}.fs]';
      carrier_cos = sqrt(2) * cos(2*pi*settings{j}.fc*t);
      carrier_sin = sqrt(2) * sin(2*pi*settings{j}.fc*t);
      %     add MFCW carriers
      signal = mfcw_addseccarriers(carrier_cos, settings{j}.mfcw_addseccarriers);
      %     transmit over channel
      signal = cell2mat( channel_main({signal}, settings{j}.channel) );
      %     modulate
      signal = signal .* (0.5+square(2*pi*settings{j}.fm*t)*settings{j}.drho/2);
      %     transmit over channel
      signal = cell2mat( channel_main({signal}, settings{j}.channel) );
      %     demodulation, AAF and downsampling (zero-delay filter)
      [b,a] = cheby2(internalsettings.aaf_ord, internalsettings.aaf_att/2, settings{j}.frs/settings{j}.fs);
      signal = complex(signal.*carrier_cos, -signal.*carrier_sin); 
      signal = filtfilt(b,a, signal);
      signal = signal(round(1 : settings{j}.fs/settings{j}.frs : end));
      %     component selection
      [c_mi, c_i, c_im] = mfcw_compsel(signal, settings{j}.mfcw_compsel);
      %     averaging
      ind = round(1+length(c_mi)*internalsettings.est_ind(1)):round(length(c_mi)*internalsettings.est_ind(2));
      avg_c_mi = mean(c_mi(ind,:));
      avg_c_im = mean(c_im(ind,:));
      %     distance estimation
      dist_hat = mfcw_calcdist(avg_c_mi, avg_c_im, settings{j}.mfcw_calcdist);
      
      % checks
      for nc = 2 : settings{j}.nc + 1
         [errors, errortext] = check_result(j, dist_hat{nc}, settings{j}.dist,...
            errors, errortext, sprintf('%iFCW',nc), 'relerr', internalsettings.tol);
      end
   end
   
   % final output for this test block
   if errors == 0
      disp('         ... passed');
   else
      disp(sprintf('         ... ERRORS: %s', errortext));
      sumoferrors = sumoferrors + 1;
   end
end


% switch warning(s) back on again
warning(oldwarn);



return
% *******************************************************************************************************
% DEBUGGING SNIPPETS

% disp(sprintf('%2i:%s', j, sprintf(' %9.2f%% ', 100*(cell2mat(dist_hat) - settings{j}.dist) / settings{j}.dist)));

% close all
% for kk = 1 : settings{j}.nc+1
%    figure; hold on;
%    plot(abs(c_mi(:,kk)), 'b');
%    plot(abs(c_im(:,kk)), 'k');
%    plot(ind, abs(c_mi(ind,kk)), 'r');
%    plot(ind, abs(c_im(ind,kk)), 'r');
%    hold off; grid on; xlim([ind(1)-ind(1)/2, ind(end)+ind(1)/2]);
% end

% nfft = 2^floor(log2(length(signal)));
% figure; pwelch(signal, rectwin(nfft),[],nfft,settings{j}.fs); grid minor;    
% figure; pwelch(signal, rectwin(nfft),[],nfft,settings{j}.frs); grid minor;
% freqz(b,a, linspace(0,settings{i}.fs/2,1e3), settings{i}.fs);

