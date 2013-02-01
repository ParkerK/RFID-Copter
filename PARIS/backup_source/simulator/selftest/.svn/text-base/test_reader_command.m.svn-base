% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: reader - command
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
% version = test_reader_command()
%    Just returns the version number (string).
% sumoferrors = test_reader_command(sumoferrors)
%    Tests the function reader_command, displays the results in the command window and returns the 
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
%
% *******************************************************************************************************


function sumoferrors = test_reader_command(sumoferrors)
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
% initialization

% output
disp('   = reader_command **');

% load expected results and settings
load results_reader_command

% load expected results and modulator settings
data = loadmat('results_reader_command.mat', '     ');
settings = data.settings;
results  = data.results;


% *******************************************************************************************************
% test: query

disp(sprintf('      - command: query (%i tests with different setups)', length(settings.query)));
% run tests and evaluate
errors = 0;
for i = 1 : length(settings.query)
   data = reader_command('query', settings.query{i});
   errors = check_result(i, data, results.query(i,:), errors, '', '', 'equal');
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERROR (%i tests ok, %i not ok)', length(settings.query)-errors, errors));
   sumoferrors = sumoferrors + 1;
end


% *******************************************************************************************************
% test: ack

disp(sprintf('      - command: ack (%i tests with different setups)', length(settings.ack)));
% run tests and evaluate
errors = 0;
for i = 1 : length(settings.ack)
   data = reader_command('ack', settings.ack{i});
   errors = check_result(i, data, results.ack(i,:), errors, '', '', 'equal');
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERROR (%i tests ok, %i not ok)', length(settings.ack)-errors, errors));
   sumoferrors = sumoferrors + 1;
end
