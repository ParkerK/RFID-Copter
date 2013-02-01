% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: channel - main file
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
% version = test_channel_main()
%    Just returns the version number (string).
% sumoferrors = test_channel_main(sumoferrors)
%    Tests the function CHANNEL_MAIN, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_channel_main(sumoferrors);
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


function sumoferrors = test_channel_main(sumoferrors)
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
% initialize

% output
disp('   = channel_main **');
disp('     checking impulse response +/-0% and noise variance +/-5% for random setups')

% load settings and expected results
data = loadmat('results_channel_main.mat', '     ');
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
   
   % display once for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
   end
   
   % reset signals
   clear('input', 'output');
          
   % create unit pulses (input signals) of different length
   for tx = 1 : settings{i}.ntx
      input{tx} = [settings{i}.channel{tx,1}.txa; zeros(settings{i}.tx_len(tx)-1, 1)];
   end
   
   % apply channel
   output = channel_main(input, settings{i}.channel);
         
   % compare to expected results
   %     reset error/errortext once for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      errors    = 0;
      errortext = '';
   end
   %     mandatory checks and conversions
   %        ... check output signal lengths => have to be identical
   if any(diff(cellfun(@length, output)) ~= 0)
      [errors, errortext] = check_result(i, true, false, errors, errortext, 'length mismatch', 'equal'); % always an error
      continue;
   end
   %        ... all RX signals have identical lengths => channel impulse responses can be represented by
   %            a matrix (simpler to handle)
   cir = cell2mat(output);
   %        ... compare length to expected
   [errors, errortext] = check_result(i, results{i}.rx_len, size(cir,1), errors, errortext, 'length', 'equal');
   %     for noise: check variance (remove expected result from signal => everything else is "noise")
   if settings{i}.switches.noise_on
      [errors, errortext] = check_result(i, var(cir-results{i}.ol), results{i}.nvar,...
         errors, errortext, 'var', 'relerr', 0.05);
   %     otherwise: check impulse response
   else
      [errors, errortext] = check_result(i, sparse(cir), results{i}.ol,...
         errors, errortext, 'impres', 'equal');
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

