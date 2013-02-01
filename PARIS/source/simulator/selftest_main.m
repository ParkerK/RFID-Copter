% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% main script for self test
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
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
%
%
% ***** Todo *****
% = bring outdated selftests up to date
%
% *******************************************************************************************************


% *******************************************************************************************************
% initialization
clear all; close all; clc;

% paths, etc
%     path to globalinit.m
%path = fullfile(fileparts(mfilename('fullpath')), '..');
fullfile(fileparts(mfilename('fullpath')), '..')
path = '/Users/johnpohollaren/Documents/Reynolds_Lab/RFID-Copter/PARIS/source/';
%     Matlab does not resolve symbolic links; if we're running on Unix, let the system resolve them
if isunix & ~ismac
  [dummy, path] = system(sprintf('readlink -f "%s"', path));
end
%     add this path and call globalinit script

addpath(path); clear('dummy', 'path');
globalinit('silent');
%     "switch to here" (also necessary for version_system('all'))
cd(fileparts(mfilename('fullpath')));

% filename of diary file
diaryfile = 'lastselftest.log';

% sum of errors
errors_reader  = 0;
errors_channel = 0;
errors_tag     = 0;
errors_ranging = 0;
errors_main    = 0;

% error structs (try/catch)
error_structs = {};


% *******************************************************************************************************
% diary and headline

% enable diary
system(sprintf('rm -f %s', diaryfile));
diary(diaryfile);

% headline
disp(sprintf('Created on %s', datestr(now, 'local')));



% *******************************************************************************************************
% version of all functions
disp(' ');
disp('*******************************************************************************************************');
disp('* Versions:');

global globalsettings;
globalsettings.logging.versions = 1;
globalsettings.logging.en_svn   = 1;
version_system('all');
globalsettings.logging.versions = 0;
globalsettings.logging.en_svn   = 0;


% *******************************************************************************************************
% reader
disp(' ');
disp('*******************************************************************************************************');
disp('* Reader:');

try errors_reader = test_reader_oscillator(errors_reader); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_reader = errors_reader + 1;
end

try errors_reader = test_reader_command(errors_reader); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_reader = errors_reader + 1;
end

try errors_reader = test_reader_modulation(errors_reader); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_reader = errors_reader + 1;
end

try errors_reader = test_reader_transmitter(errors_reader); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_reader = errors_reader + 1;
end

try errors_reader = test_reader_receiver(errors_reader); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_reader = errors_reader + 1;
end

try errors_reader = test_reader_demodulation(errors_reader); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_reader = errors_reader + 1;
end

try errors_reader = test_reader_analogpathest(errors_reader); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_reader = errors_reader + 1;
end


% *******************************************************************************************************
% channel
disp(' ');
disp('*******************************************************************************************************');
disp('* Channel:');

% try errors_channel = test_channel_newpos(errors_channel); catch me
%    disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
%    error_structs = [error_structs, me];
%    errors_channel = errors_channel + 1;
% end

try errors_channel = test_channel_large(errors_channel); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_channel = errors_channel + 1;
end

try errors_channel = test_channel_directivity(errors_channel); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_channel = errors_channel + 1;
end

try errors_channel = test_channel_surface(errors_channel); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_channel = errors_channel + 1;
end

try errors_channel = test_channel_small(errors_channel); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_channel = errors_channel + 1;
end

try errors_channel = test_channel_noise(errors_channel); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_channel = errors_channel + 1;
end

% try errors_channel = test_channel_main(errors_channel); catch me
%    disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
%    error_structs = [error_structs, me];
%    errors_channel = errors_channel + 1;
% end

disp(' ... temporarily deactivated tests (outdated): channel_newpos, channel_main');


% *******************************************************************************************************
% tag
disp(' ');
disp('*******************************************************************************************************');
disp('* Tag:');

try errors_tag = test_tag_clock(errors_tag); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_tag = errors_tag + 1;
end

try errors_tag = test_tag_encoding(errors_tag); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_tag = errors_tag + 1;
end

try errors_tag = test_tag_decoding(errors_tag); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_tag = errors_tag + 1;
end

try errors_tag = test_tag_modulation(errors_tag); catch me
   disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
   error_structs = [error_structs, me];
   errors_tag = errors_tag + 1;
end


% *******************************************************************************************************
% ranging
disp(' ');
disp('*******************************************************************************************************');
disp('* Ranging:');

disp(' ... this selftest is outdated and has been deactivated temporarily');

% try errors_ranging = test_mfcw(errors_ranging); catch me
%    disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
%    error_structs = [error_structs, me];
%    errors_ranging = errors_ranging + 1;
% end


% *******************************************************************************************************
% main
disp(' ');
disp('*******************************************************************************************************');
disp('* Main:');

disp(' ... this selftest is outdated and has been deactivated temporarily');

% try errors_main = test_main_reader_tag(errors_main); catch me
%    disp(sprintf('     ... FAILED WITH AN ERROR (see %d. below)', length(error_structs)+1));
%    error_structs = [error_structs, me];
%    errors_main = errors_main + 1;
% end


% *******************************************************************************************************
% print caught errors
if ~isempty(error_structs)
   disp(' ');
   disp('*******************************************************************************************************');
   disp('* THERE WERE SOME EXCEPTIONS:');
   for i = 1 : length(error_structs)
      disp(sprintf('\n%2d. %s line %d: "%s"\n    %s', i, error_structs(i).stack(1).name,...
         error_structs(i).stack(1).line, error_structs(i).identifier,...
         strrep(error_structs(i).message, sprintf('\n'), sprintf('\n    ')) ));
   end
end


% *******************************************************************************************************
% result
disp(' ');
disp('*******************************************************************************************************');
disp('* RESULT:');
disp(sprintf('Reader:  %3i error(s)\nChannel: %3i error(s)\nTag:     %3i error(s)\nRanging: %3i error(s)\nMain:    %3i error(s)',...
   errors_reader, errors_channel, errors_tag, errors_ranging, errors_main));
if errors_reader + errors_channel + errors_tag + errors_ranging + errors_main > 0
   disp(sprintf('\n   ... note that false positive errors are possible due to the statistical nature of most tests.'));
end
disp(sprintf('\nLegend:'));
disp('**  statistical test (not all possible setup combinations covered)');
disp('*** low test coverage (not all parameters tested)');


% *******************************************************************************************************
% cleanup

diary off;

