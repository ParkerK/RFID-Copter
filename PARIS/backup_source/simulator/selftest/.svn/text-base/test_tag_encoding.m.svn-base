% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: tag - encoding (including CRC-16 calculation)
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
% version = test_tag_encoding()
%    Just returns the version number (string).
% sumoferrors = test_tag_encoding(sumoferrors)
%    Tests the function tag_encoding, displays the results in the command window and returns the overall
%    amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_tag_encoding(sumoferrors);
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



function sumoferrors = test_tag_encoding(sumoferrors)
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
disp('   = tag_encoding');

% load expected results
data = loadmat('results_tag_encoding', '     ');
results  = data.results;


% *******************************************************************************************************
% test: FM0 encoding, no pilot

disp('      - FM0, TRext=0 (no pilot tone)');
% settings
settings.tag.m     = 1;
settings.tag.trext = 0;
% initialize test result variables
errors = 0;
errortext = '';
% run tests and evaluate
for i = 0 : 3 % for all sequences
   encoded = tag_encoding(bitget(i, 2:-1:1), settings.tag);
   [errors, errortext] = check_result(i, encoded, results.fm0.trext0(:,i+1), errors, errortext, dec2bin(i,2));
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS in sequence(s)%s', errortext));
   sumoferrors = sumoferrors + 1;
end



% *******************************************************************************************************
% test: FM0 encoding, pilot

disp('      - FM0, TRext=1 (pilot tone)');
% settings
settings.tag.m     = 1;
settings.tag.trext = 1;
% initialize test result variables
errors = 0;
errortext = '';
% run tests and evaluate
for i = 0 : 3 % for all sequences
   encoded = tag_encoding(bitget(i, 2:-1:1), settings.tag);
   [errors, errortext] = check_result(i, encoded, results.fm0.trext1(:,i+1), errors, errortext, dec2bin(i,2));
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS in sequence(s)%s', errortext));
   sumoferrors = sumoferrors + 1;
end


% *******************************************************************************************************
% test: Miller encoding, M=2, no pilot

disp('      - Miller, M=2, TRext=0 (no pilot tone)');
% settings
settings.tag.m        = 2;
settings.tag.trext    = 0;
% initialize test result variables
errors = 0;
errortext = '';
% run tests and evaluate
for i = 0 : 7 % for all sequences
   encoded = tag_encoding(bitget(i, 3:-1:1), settings.tag);
   [errors, errortext] = check_result(i, encoded, results.miller2.trext0(:,i+1), errors, errortext, dec2bin(i,3));
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS in sequence(s)%s', errortext));
   sumoferrors = sumoferrors + 1;
end


% *******************************************************************************************************
% test: Miller encoding, M=2, pilot

disp('      - Miller, M=2, TRext=1 (pilot tone)');
% settings
settings.tag.m        = 2;
settings.tag.trext    = 1;
% initialize test result variables
errors = 0;
errortext = '';
% run tests and evaluate
for i = 0 : 7 % for all sequences
   encoded = tag_encoding(bitget(i, 3:-1:1), settings.tag);
   [errors, errortext] = check_result(i, encoded, results.miller2.trext1(:,i+1), errors, errortext, dec2bin(i,3));
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS in sequence(s)%s', errortext));
   sumoferrors = sumoferrors + 1;
end


% *******************************************************************************************************
% test: Miller encoding, M=4, no pilot

disp('      - Miller, M=4, TRext=0 (no pilot tone)');
% settings
settings.tag.m        = 4;
settings.tag.trext    = 0;
% initialize test result variables
errors = 0;
errortext = '';
% run tests and evaluate
for i = 0 : 7 % for all sequences
   encoded = tag_encoding(bitget(i, 3:-1:1), settings.tag);
   [errors, errortext] = check_result(i, encoded, results.miller4.trext0(:,i+1), errors, errortext, dec2bin(i,3));
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS in sequence(s)%s', errortext));
   sumoferrors = sumoferrors + 1;
end


% *******************************************************************************************************
% test: Miller encoding, M=4, pilot

disp('      - Miller, M=4, TRext=1 (pilot tone)');
% settings
settings.tag.m        = 4;
settings.tag.trext    = 1;
% initialize test result variables
errors = 0;
errortext = '';
% run tests and evaluate
for i = 0 : 7 % for all sequences
   encoded = tag_encoding(bitget(i, 3:-1:1), settings.tag);
   [errors, errortext] = check_result(i, encoded, results.miller4.trext1(:,i+1), errors, errortext, dec2bin(i,3));
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS in sequence(s)%s', errortext));
   sumoferrors = sumoferrors + 1;
end


% *******************************************************************************************************
% test: Miller encoding, M=8, no pilot

disp('      - Miller, M=8, TRext=0 (no pilot tone)');
% settings
settings.tag.m        = 8;
settings.tag.trext    = 0;
% initialize test result variables
errors = 0;
errortext = '';
% run tests and evaluate
for i = 0 : 7 % for all sequences
   encoded = tag_encoding(bitget(i, 3:-1:1), settings.tag);
   [errors, errortext] = check_result(i, encoded, results.miller8.trext0(:,i+1), errors, errortext, dec2bin(i,3));
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS in sequence(s)%s', errortext));
   sumoferrors = sumoferrors + 1;
end


% *******************************************************************************************************
% test: Miller encoding, M=8, pilot

disp('      - Miller, M=8, TRext=1 (pilot tone)');
% settings
settings.tag.m        = 8;
settings.tag.trext    = 1;
% initialize test result variables
errors = 0;
errortext = '';
% run tests and evaluate
for i = 0 : 7 % for all sequences
   encoded = tag_encoding(bitget(i, 3:-1:1), settings.tag);
   [errors, errortext] = check_result(i, encoded, results.miller8.trext1(:,i+1), errors, errortext, dec2bin(i,3));
end
% output test result(s)
if errors == 0
   disp('         ... passed');
else
   disp(sprintf('         ... ERRORS in sequence(s)%s', errortext));
   sumoferrors = sumoferrors + 1;
end

