% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: tag - clock (sample time vector generator)
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
% version = test_tag_clock()
%    Just returns the version number (string).
% sumoferrors = test_tag_clock(sumoferrors)
%    Tests the function tag_clock, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_tag_clock(sumoferrors);
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


function sumoferrors = test_tag_clock(sumoferrors)
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
disp('   = tag_clock **');

% load expected results
data = loadmat('results_tag_clock', '     ');
settings  = data.settings;


% *******************************************************************************************************
% test: mode = 'fs'

disp(sprintf('      - %i different setups in mode "fs"', length(settings.mode_fs)));
disp('        checking length, phi0, fcenter (+/-0.1%), fsigma (+/-5%)'); 

% reset error counters
errors.length  = 0;
errors.phi0    = 0;
errors.fcenter = 0;
errors.fsigma  = 0;

% run tests and evaluate
for i = 1 : length(settings.mode_fs)
   % generate vector
   testsetup = settings.mode_fs(i);
   rand('twister', i);  % initialize RNG for reproducible results
   clock = tag_clock(testsetup);
    
   % obtain parameters
   %     phi0
   t0   = clock(1) - 1; 
   phi0 = t0 * testsetup.fcenter / testsetup.fs * 360;
   %     fcenter, fsigma: fit linear model to derive clean clock
   n = [1:1:length(clock)]';
   lmodel = polyfit(n, clock-t0, 1);
   clean  = polyval(lmodel, n);
   %     fcenter
   fcenter = testsetup.fs / lmodel(1);
   %     fsigma (use settings.fcenter to minimize variance)
   fsigma = std(clock-clean) * testsetup.fcenter / testsetup.fs * 100;
   
   % check parameters
   errors.length    = max(clock) > testsetup.length;
   errors.phi0      = check_result(i, phi0,    testsetup.phi0_res, errors.phi0,    '', '', 'range', 180*testsetup.fcenter/testsetup.fs);
   errors.fcenter   = check_result(i, fcenter, testsetup.fcenter,  errors.fcenter, '', '', 'relerr', 0.001);
   errors.fsigma    = check_result(i, fsigma,  testsetup.fsigma,   errors.fsigma,  '', '', 'relerr', 0.05); % mse=0.05 if testsetup.fsigma==0 
end

% output test result(s)
if errors.length + errors.phi0 + errors.fcenter + errors.fsigma == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS: %i length, %i phi0, %i fcenter, %i fsigma',...
            errors.length, errors.phi0, errors.fcenter, errors.fsigma));
   sumoferrors = sumoferrors + 1;
end


% *******************************************************************************************************
% test: mode = 'fclk'

disp(sprintf('      - %i different setups in mode "fclk"', length(settings.mode_fs)));
disp('        checking length, phi0, fcenter (+/-0.1%), fsigma (+/-5%)'); 

% reset error counters
errors.length  = 0;
errors.phi0    = 0;
errors.fcenter = 0;
errors.fsigma  = 0;

% run tests and evaluate
for i = 1 : length(settings.mode_fclk)
   % generate vector
   testsetup = settings.mode_fclk(i);
   rand('twister', 12345+i);  % initialize RNG for reproducible results
   clock = tag_clock(testsetup);
    
   % obtain parameters
   %     phi0
   t0   = clock(1) - 1; 
   phi0 = t0 * testsetup.fcenter / testsetup.fs * 360;
   %     fcenter, fsigma: fit linear model to derive clean clock
   n = [1:1:length(clock)]';
   lmodel = polyfit(n, clock-t0, 1);
   clean  = polyval(lmodel, n);
   %     fcenter
   fcenter = testsetup.fs / lmodel(1);
   %     fsigma (use settings.fcenter to minimize variance)
   fsigma = std(clock-clean) * testsetup.fcenter / testsetup.fs * 100;
   
   % check parameters
   errors.length  = length(clock) ~= testsetup.length;
   errors.phi0    = check_result(i, phi0,    testsetup.phi0_res, errors.phi0,    '', '', 'range', 180*testsetup.fcenter/testsetup.fs);
   errors.fcenter = check_result(i, fcenter, testsetup.fcenter,  errors.fcenter, '', '', 'relerr', 0.001);
   errors.fsigma  = check_result(i, fsigma,  testsetup.fsigma,   errors.fsigma,  '', '', 'relerr', 0.05); % mse=0.05 if testsetup.fsigma==0 
end

% output test result(s)
if errors.length + errors.phi0 + errors.fcenter + errors.fsigma == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS: %i length, %i phi0, %i fcenter, %i fsigma',...
            errors.length, errors.phi0, errors.fcenter, errors.fsigma));
   sumoferrors = sumoferrors + 1;
end
