% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: reader - modulation
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
% version = test_reader_modulation()
%    Just returns the version number (string).
% sumoferrors = test_reader_modulation(sumoferrors)
%    Tests the function reader_modulation, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_reader_modulation(sumoferrors);
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
% = add SSB-ASK after overshoot problem is fixed
% = time-variant checks (moddepth, ...)
%
% *******************************************************************************************************


function sumoferrors = test_reader_modulation(sumoferrors)
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
%     avg times/amplitudes
internalsettings.tol_avg = 0.01; % relative error
internalsettings.tol_t0  = 0.05; % relative error (t0)
%     standard deviations
internalsettings.tol_std = 0.05; % times average value


% *******************************************************************************************************
% initialization

% output
disp('   = reader_modulation **');
disp(sprintf('     checking avg times/amplitudes +/-%g%% (+/-%g%% for t0) and stddev +/-%g%% of avg',...
   internalsettings.tol_avg*100, internalsettings.tol_t0*100, internalsettings.tol_std*100));
disp('     WARNING: no ripple or time-variant envelope checks (e.g. for moddepth) done');

% load expected results and modulator settings
data = loadmat('results_reader_modulation.mat', '     ');
settings     = data.settings;
partitioning = data.partitioning;
oscsettings  = data.oscsettings;
results      = data.results;

% create carrier
bigcarrier = reader_oscillator(oscsettings);


% *******************************************************************************************************
% run tests and evaluate

for i = 1 : length(results)
   
   % partitioning and output
   %     find index for partition struct
   for j = 1 : length(partitioning.indices)
      if i <= partitioning.indices(j)
         index = j;
         break;
      end
   end
   
   % display once for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
   end
   
   % initialize RNGs
   rand('twister', i);
   randn('state', i);
   
   % determine needed carrier length, complete and sanitize settings 
   settings{i} = reader_modulation(settings{i}.data, settings{i});
   
   % compare created settings{i} to results{i}
   set = struct2cell(settings{i});
   res = struct2cell(results{i});
   if sum(cellfun(@isequal, set, res)) ~= length(set)
      disp('         ... ERROR: created settings not equal to expected settings; skipping modulation');
      sumoferrors = sumoferrors + 1;
      continue;
   end
   
   % modulate
   modcarrier = reader_modulation(bigcarrier(1:settings{i}.length_s), settings{i}.data, settings{i});
   
   % measure timings
   timings = linktiming_reader(modcarrier);
     
   % compare to expected results
   %     reset error/errortext once for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      errors    = 0;
      errortext = '';
   end
   %     check timing values
   [errors, errortext] = check_result(i, timings.moddepth(1), settings{i}.moddepth,...
      errors, errortext, 'moddepth',  'relerr', internalsettings.tol_avg);
   [errors, errortext] = check_result(i, timings.t0_s,        settings{i}.t0_s,...
      errors, errortext, 't0',        'relerr', internalsettings.tol_t0);
   [errors, errortext] = check_result(i, timings.delimiter_s, settings{i}.delimiter_s,...
      errors, errortext, 'delimiter', 'relerr', internalsettings.tol_avg);
   [errors, errortext] = check_result(i, timings.tari_s,      settings{i}.tari_s,...
      errors, errortext, 'tari',      'relerr', internalsettings.tol_avg);
   [errors, errortext] = check_result(i, timings.rtcal_s,     settings{i}.rtcal_s,...
      errors, errortext, 'rtcal',     'relerr', internalsettings.tol_avg);
   if strcmpi(settings{i}.leadin, 'preamble') % no trcal for framesync
      [errors, errortext] = check_result(i, timings.trcal_s,  settings{i}.trcal_s,...
         errors, errortext, 'trcal',     'relerr', internalsettings.tol_avg);
   end
   [errors, errortext] = check_result(i, timings.x_s(1),  settings{i}.x_s,...
      errors, errortext, 'x',     'relerr', internalsettings.tol_avg);
   [errors, errortext] = check_result(i, timings.pw_s(1), settings{i}.pw_s,...
      errors, errortext, 'pw',    'relerr', internalsettings.tol_avg);
   [errors, errortext] = check_result(i, timings.tr_s(1), settings{i}.trf_s,...
      errors, errortext, 'tr',    'relerr', internalsettings.tol_avg);
   [errors, errortext] = check_result(i, timings.tf_s(1), settings{i}.trf_s,...
      errors, errortext, 'tf',    'relerr', internalsettings.tol_avg);
   [errors, errortext] = check_result(i, timings.data0_s(1), settings{i}.tari_s,...
      errors, errortext, 'data0', 'relerr', internalsettings.tol_avg);
   [errors, errortext] = check_result(i, timings.data1_s(1), settings{i}.tari_s+settings{i}.x_s,...
      errors, errortext, 'data1', 'relerr', internalsettings.tol_avg);
   %     check standard deviations
   [errors, errortext] = check_result(i, timings.moddepth(2), 0,...
      errors, errortext, 'stddev moddepth', 'range', internalsettings.tol_std*timings.moddepth(1));
   [errors, errortext] = check_result(i, timings.x_s(2),      0,...
      errors, errortext, 'stddev x',        'range', internalsettings.tol_std*timings.x_s(1));
   [errors, errortext] = check_result(i, timings.pw_s(2),     0,...
      errors, errortext, 'stddev pw',       'range', internalsettings.tol_std*timings.pw_s(1));
   [errors, errortext] = check_result(i, timings.tr_s(2),     0,...
      errors, errortext, 'stddev tr',       'range', internalsettings.tol_std*timings.tr_s(1));
   [errors, errortext] = check_result(i, timings.tf_s(2),     0,...
      errors, errortext, 'stddev tf',       'range', internalsettings.tol_std*timings.tf_s(1));
   [errors, errortext] = check_result(i, timings.data0_s(2),  0,...
      errors, errortext, 'stddev data0',    'range', internalsettings.tol_std*timings.data0_s(1));
   [errors, errortext] = check_result(i, timings.data1_s(2),  0,...
      errors, errortext, 'stddev data1',    'range', internalsettings.tol_std*timings.data1_s(1));
   
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

