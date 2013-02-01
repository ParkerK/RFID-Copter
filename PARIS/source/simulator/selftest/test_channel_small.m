% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: channel - small scale
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
% version = test_channel_small()
%    Just returns the version number (string).
% sumoferrors = test_channel_small(sumoferrors)
%    Tests the function channel_small, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_channel_small(sumoferrors);
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
% ? narrowband correlations
%
% *******************************************************************************************************


function sumoferrors = test_channel_small(sumoferrors)
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

% relative error tolerances for the individual tests
internalsettings.tol_k    = 1e-2; % rel. error K-factor
internalsettings.tol_pwr  = 1e-2; % rel. error LOS power
internalsettings.tol_trms = 5e-2; % rel. error RMS delay spread
internalsettings.tol_k_stat    = internalsettings.tol_k/1e6; % tolerance in rel. error K-factor in statistics (numerics)
internalsettings.tol_trms_stat = internalsettings.tol_trms/1e6; % tolerance in rel. error RMS delay spread in statistics

% size for randomization check (channel_small sets seed for rand and randn; make sure they are randomized again)
internalsettings.len_randchk = 10; % samples


% *******************************************************************************************************
% initialize

% output
disp('   = channel_small **');
disp(sprintf('     checking #paths,iterations<=max; avg./est. PDP: trms +/-%g%%, LOS pwr +/-%g%%, K +/-%g%%',...
   internalsettings.tol_pwr*100, internalsettings.tol_trms*100, internalsettings.tol_k*100));

% load settings and expected results
data = loadmat('results_channel_small.mat', '     ');
settings     = data.settings;
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
   
   % once for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      disp(sprintf('      - %i setup(s) "%s", %g ensemble(s) per setup',...
         partitioning.runs(index), partitioning.names{index}, settings{i}.ensembles));
      if ~strcmpi(partitioning.names{index}, 'off')
         if any(strcmpi(partitioning.names{index}, {'deterministic NLOS + seed', 'random NLOS + seed'}))
            disp('        checking rand/randn after channel_small (randomization) and reproducability (identical channels)');
         else
            disp('        checking rand/randn after channel_small (randomization)');
         end
      end
      % initialize randomization check matrices
      seedtest.rand    = [];
      seedtest.randn   = [];
      seedtest.gain    = {};
      seedtest.delay_s = {};
   end
   
   % call channel_small (plus ensemble-averaging)
   [delay_s, gain, statistics] = channel_small(settings{i});
   %     average power => overall gain
   gain = sqrt( sum(gain.^2, 2) / settings{i}.ensembles );
   %     average delay vector
   delay_s = mean(delay_s, 2);
      
   % random sequence tests (for all settings)
   seedtest.rand  = [seedtest.rand ; rand( 1,internalsettings.len_randchk)];
   seedtest.randn = [seedtest.randn; randn(1,internalsettings.len_randchk)];
   
   % reproducability tests: record gain/delay (lengths might not match in case of an error => cells)
   if ~isnan(settings{i}.seed)
      seedtest.gain    = [seedtest.gain,    gain];
      seedtest.delay_s = [seedtest.delay_s, delay_s];
   end
   
   % K: [dB] -> linear
   settings{i}.k_lin = 10^(settings{i}.k/10);
   statistics.av_k_lin = 10^(statistics.av_k/10);
        
   % compare to expected results
   %     reset error/errortext once for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      errors    = 0;
      errortext = '';
   end   
   
   %     if small-scale model was switched on
   if settings{i}.on
      % maximum lengths / etc
      [errors, errortext] = check_result(i, statistics.iter, settings{i}.maxiter,...
         errors, errortext, 'maxiter', 'range', [-Inf, 0]); % smaller or equal
      %     forced # of taps will overwrite maxrays, but not the removal of identical entries
      %     => don't check if a resolution warning has been issued
      if isfield(settings{i}, 'nrays_f') && ~statistics.warn_res
         [errors, errortext] = check_result(i, length(delay_s), settings{i}.nrays_f,...
            errors, errortext, 'forced # rays', 'equal');
      elseif ~isfield(settings{i}, 'nrays_f')
         [errors, errortext] = check_result(i, length(delay_s), settings{i}.maxrays,...
            errors, errortext, 'maxrays', 'range', [-Inf, 0]); % smaller or equal
      end
      
      % comparison settings to reported statistics
      % ... only if iteration was complete, otherwise the mismatch can be arbitrary
      if ~(statistics.warn_res || statistics.warn_bw || statistics.warn_ktrms ||...
           statistics.warn_maxiter || statistics.warn_maxrays)
         [errors, errortext] = check_result(i, statistics.av_trms, settings{i}.trms,...
            errors, errortext, 'av_trms', 'relerr', settings{i}.eps_trms*(1+internalsettings.tol_trms_stat));
         [errors, errortext] = check_result(i, statistics.av_k_lin, settings{i}.k_lin,...
            errors, errortext, 'av_K', 'relerr', settings{i}.eps_k*(1+internalsettings.tol_k_stat));
      end
      
      % estimates to reported statistics: timings
      %     trms
      [errors, errortext] = check_result(i, sqrt(var(delay_s/settings{i}.fs, gain.^2)), statistics.av_trms,...
         errors, errortext, 'trms', 'relerr', internalsettings.tol_trms);
      % estimates to reported statistics: gain/power
      %     Ricean K-Factor
      [errors, errortext] = check_result(i, gain(1)^2/sum(gain(2:end).^2), statistics.av_k_lin,...
         errors, errortext, 'K/eps', 'relerr', internalsettings.tol_k);
      %     LOS power
      [errors, errortext] = check_result(i, gain(1)^2, statistics.av_k_lin/(statistics.av_k_lin+1),...
         errors, errortext, 'LOS', 'relerr', internalsettings.tol_pwr);
      
      % if this is the last test in block
      if i == partitioning.indices(index)
         % check if all random sequences after the call of channel_small are different
         [errors, errortext] = check_result(i, all(sum(diff(seedtest.rand)==0,2)), 0,...
            errors, errortext, 'randomization after return (rand)',  'equal');
         [errors, errortext] = check_result(i, all(sum(diff(seedtest.randn)==0,2)), 0,...
            errors, errortext, 'randomization after return (randn)', 'equal');
         % if this was a seed test (reproducible results): check if all gains/delays are identical
         if ~isnan(settings{i}.seed)
            % check lengths first
            [errors, errortext] = check_result(i, all(diff(cellfun(@(x)length(x), seedtest.gain))==0), 1,...
               errors, errortext, 'repeatability (length gain)',  'equal');
            [errors, errortext] = check_result(i, all(diff(cellfun(@(x)length(x), seedtest.delay_s))==0), 1,...
               errors, errortext, 'repeatability (length delay_s)', 'equal');
            % if lengths match: are all the contents identical 
            if all(diff(cellfun(@(x)length(x), seedtest.gain))==0) && all(diff(cellfun(@(x)length(x), seedtest.delay_s))==0)
               [errors, errortext] = check_result(i, sum(sum(diff(cell2mat(seedtest.gain),1,2))), 0,...
                  errors, errortext, 'repeatability (gain)',  'equal');
               [errors, errortext] = check_result(i, sum(sum(diff(cell2mat(seedtest.delay_s),1,2))), 0,...
                  errors, errortext, 'repeatability (delay_s)', 'equal');
            end   
         end
      end
   
   %     if small-scale model was switched off
   else
      [errors, errortext] = check_result(i, [delay_s, gain], [0, 1], errors, errortext, 'mode:off', 'equal');
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


%          [errors, errortext] = check_result(i, abs(1-statistics.av_k_lin/settings{i}.k_lin), settings{i}.eps_k,...
%             errors, errortext, 'av_K', 'range', [-Inf, 0]); % smaller or equal

%          [errors, errortext] = check_result(i, abs(1-statistics.av_k_lin/settings{i}.k_lin), settings{i}.eps_k,...
%             errors, errortext, 'av_K', 'range', [-Inf, internalsettings.tol_k(index)*settings{i}.eps_k]);
