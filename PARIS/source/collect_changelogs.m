% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% Collect changelogs for all files that are managed by svn and have been changed in some way and writes
% them to file (using diary). This file can then be used with the commit command.
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
% - more performant implementation using cellfun
%
% *******************************************************************************************************


% *******************************************************************************************************
% initialization and settings
clear all; close all; clc;
version = 'beta 3.0';

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% initialize simulator/emulator version control (also used by this script)
addpath(fullfile(fileparts(mfilename('fullpath'))));
globalinit('init');

% for search
settings.rootpath = fileparts(mfilename('fullpath'));
settings.goodstatus  = '%[ACMR]'; % status characters to search for (formated for textread)

% to select changelog out of header
settings.startmarker = '% ***** Changelog *****';
settings.midmarker   = '% ***** Todo *****';
settings.endmarker   = '% *******************************************************************************************************';

% other
settings.filename = 'changes.log'; % only relativ paths!
settings.commitcommand = 'svn commit --file %s';


% *******************************************************************************************************
% check if svn is available (otherwise the rest does not make sense) and if so switch directory

if ~svn_isavailable()
   disp('Sorry, svn not available.');
   return
end

% switch directory to rootpath
cd(settings.rootpath);


% *******************************************************************************************************
% get list of all files that have status not equal to (recursively)

% get a list of all files recursively (only .m-files have a changelog)
files = rdir(fullfile(settings.rootpath, '**/*.m'));
files = {files(:).name}'; % we're only interested in the names

% remove all files that have not status = goodstatus
goodfiles = cell(0,1);
badfiles  = cell(0,1); % files not under version control
for i = 1 : length(files)
   % get status from svn
   status = svn_status(files{i});
   % error in case a file is not under version control
   if status.svn_notcontrolled
      badfiles = [badfiles, files{i}];
   end
   % check status
   charcheck = textscan(status.svn_status, settings.goodstatus);
   % if charcheck{:} is not empty: status contaings characters from settings.goodstatus => add files
   if ~isempty(charcheck{:})
      goodfiles = [goodfiles, files{i}];
   end
end

% no files with changes
if isempty(goodfiles)
   disp('No files with changes found ...');
   return
end


% *******************************************************************************************************
% initialize diary

% start diary
system(sprintf('rm -f %s', settings.filename));
diary(settings.filename);

% add own version number
disp(sprintf('(created by Matlab script collect_changelogs, version %s on %s local time)', version, datestr(now, 'local')));


% *******************************************************************************************************
% modified files

disp(sprintf('\n\n*******************************************************************************************************'));
disp('*** CHANGELOG OF MODIFIED FILES');

for i = 1 : length(goodfiles)
   % read file line per line, 
   text = textread(goodfiles{i}, '%s', 'delimiter', '\n', 'whitespace', '');
   
   % remove leading and trailing whitespaces
   text = cellfun(@strtrim, text, 'UniformOutput', false);
   
   % find index with settings.startmarker (start of changelog)
   startline = 1;
   for line = 1 : length(text)
      if strcmp(text{line}, settings.startmarker) == 1
         startline = line;
         break;
      end
   end
   
   % find index with settings.midmarker (start of todo list)
   midline = 1;
   for line = startline : length(text)
      if strcmp(text{line}, settings.midmarker) == 1
         midline = line;
         break;
      end
   end
       
   % diplay the changelog (if present)
   if (startline > 1) && (startline+1 <= midline-2)
      disp(sprintf('\n*** %s', goodfiles{i}))
      for line = startline+1 : midline-2
         disp(text{line}(3:end))
      end
      disp(sprintf(' '));
   end
end


% *******************************************************************************************************
% todo lists

disp(sprintf('\n\n\n*******************************************************************************************************'));
disp('*** TODO LISTS OF ALL FILES');

todolists = 0;
for i = 1 : length(files)
   % read file line per line, 
   text = textread(files{i}, '%s', 'delimiter', '\n', 'whitespace', '');
   
   % remove leading and trailing whitespaces
   text = cellfun(@strtrim, text, 'UniformOutput', false);
     
   % find index with settings.midmarker (start of todo list)
   midline = 1;
   for line = 1 : length(text)
      if strcmp(text{line}, settings.midmarker) == 1
         midline = line;
         break;
      end
   end
   
   % find index with settings.endmarker (end of header)
   endline = 1;
   for line = midline : length(text)
      if strcmp(text{line}, settings.endmarker) == 1
         endline = line;
         break;
      end
   end 
         
   % display the todo list (if present)
   if (midline > 1) && (midline+1 <= endline-2)
      todolists = todolists + 1;
      disp(sprintf('\n*** %s', files{i}))
      for line = midline+1 : endline-2
         disp(text{line}(3:end))
      end
      disp(sprintf(' '));
   end
end


% *******************************************************************************************************
% files not under version control

if ~isempty(badfiles)
   disp(sprintf('\n\n\n*******************************************************************************************************'));
   disp('*** FILES NOT UNDER VERSION CONTROL');

   for i = 1 : length(badfiles)
      disp(badfiles{i})
   end
   disp(sprintf('\n\n'));
end


% ******************************************************************************************************
% other

% stop diary
diary off

% output the svn commit command (for convenience)
disp(sprintf('\n\n%i of %i files added, %i file(s) not under version control, %i files with todo lists',...
   length(goodfiles), length(files), length(badfiles), todolists));
disp(sprintf('\nJust for convenience ... here is the commit command:'));
disp(sprintf(settings.commitcommand, fullfile(pwd,settings.filename)));

