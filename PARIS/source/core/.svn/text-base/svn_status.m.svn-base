% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% call svn to get the status of a function
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
% status = svn_status(file)
%    Tries to get the svn-status of FILE and returns the information in a struct.
% status = svn_status(file, status)
%    Tries to get the svn-status of FILE and adds the information to struct STATUS
%    (does only change its own fields of the struct).
%
%
% ***** Global Variables *****
% globalsettings
%    .core.svn_head       status compared to HEAD (true) or BASE (false)
%    .core.svn_en_cache   enable SVN cache; SVN will only be called once per file if enabled
% svn_cache        cache for SVN status replies (has to be enabled)
%    .extension1      all files with extension1 (e.g. m)
%       .xyz            SVN information of file xyz.extension1 (cf. input parameter STATUS)
%    .extension2      all files with extension2 (e.g. mat)
%       .xyz            SVN information of file xyz.extension2 (cf. input parameter STATUS)
%
%
% ***** Interface definition *****
% function status = svn_status(file, status)
%    file                    filename including full path of file to check
%    status                  (optional) struct containing at least the fields stated below
%       .svn_isok            svn available and working (see SVN_ISAVAILABLE)
%       .svn_version         version of svn (see SVN_ISAVAILABLE)  
%       .svn_status          status returned by svn (first 7 characters)
%       .svn_rev             current working revision of the file
%       .oodmessage          a string like 'out of date' if file is out-of-date; empty string if not
%       .svn_notcontrolled   1 if file is NOT under version control, 0 otherwise
%       .notfound            true if file has not been found by system ("does not exist")
%
%    status                  struct containing the information obtained from svn (cf. input parameter)
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
% - recheck status from time to time; issue a warning if something has changed
%
% *******************************************************************************************************

function status = svn_status(file, status)


% *******************************************************************************************************
% settings

% global settings
global globalsettings;
global svn_cache;

% svn command to get status ... compare to BASE (default) or HEAD ?
%     The --quiet command will prevent svn_status to find files that are not under version control. In 
%     such cases status.svn_notcontrolled will be set to true.
if isstruct(globalsettings) && isstruct(globalsettings.core) &&...
      isfield(globalsettings.core, 'svn_head') && globalsettings.core.svn_head
   svn_command = 'svn status --show-updates --verbose --quiet --non-interactive'; % HEAD
else
   svn_command = 'svn status --verbose --quiet --non-interactive'; % BASE
end

% required fields in status and their init values (if status not provided)
status_requiredfields = {'svn_isok', 'svn_version', 'svn_status', 'svn_rev', 'svn_oodmessage', 'svn_notcontrolled', 'notfound'};
status_init   = {0, '', ' ', '-', '', 0, 0};

% valid status characters (formated for textscan)
validstatuschars = '%[ ACDIMRX?!~L+SKOTB]';


% *******************************************************************************************************
% input parameters

% number of input arguments
if nargin == 1
   status = cell2struct(status_init, status_requiredfields, 2);
end
if nargin < 1
   error('Not enough input arguments.')
end

% check if status contains all required fields
if ~all(isfield(status, status_requiredfields))
   error('Required fields missing in input struct "status"');
end

% (re-)initialize required fields
for i = 1 : length(status_requiredfields)
   status.(status_requiredfields{i}) = status_init{i};
end


% *******************************************************************************************************
% cache lookup

% split file in name/path/extension
[pathstr, filename, fileext] = fileparts(file);
%      remove dot from extension
fileext(1) = [];

% does the cache exist yet?
if globalsettings.core.svn_en_cache && ~isstruct(svn_cache)
   svn_cache = struct();
end

% do we have this file in cache => nothing to do for svn (?)
if globalsettings.core.svn_en_cache && isfield(svn_cache, fileext) && isfield(svn_cache.(fileext), filename)
   for i = 1 : length(status_requiredfields)
      status.(status_requiredfields{i}) = svn_cache.(fileext).(filename).(status_requiredfields{i});
   end
   return;
end


% *******************************************************************************************************
% get information from svn

% return if svn is not available (do not change anything)
[status.svn_isok, status.svn_version] = svn_isavailable();
if ~status.svn_isok
   return
end

% get file including full path
file_fullpath = which(file);
%      file does not exist => can't be controlled by SVN
if isempty(file_fullpath)
   status.notfound = 1;
   return;
end

% call svn
[dummy, reply] = system(sprintf('%s %s', svn_command, file_fullpath));

% if reply is empty: we did not search in any svn working directory or file is not under version control
% (let the calling function deal with this problem)
if isempty(reply)
   status.svn_notcontrolled = 1;
     
% otherwise: get relevant data
else
   status.svn_notcontrolled = 0;
   % split up reply in lines
   reply = textscan(reply, '%s', 'delimiter', '\n', 'Whitespace', '');
   reply = reply{:};

   % if svn found this file
   % ... if file already in svn: last line is "Status against revision:"
   % ... if file will be added: only one line, but first char is 'A'
   if ~strcmp(reply{1}(1:7), 'Status against revision')
      % characters 1 to 7 are status information
      status.svn_status = reply{1}(1:7); % complete status information
      % character 8 is out-of-date information ('*': out-of-date)
      if reply{1}(8) == '*'
         status.svn_oodmessage = ', out-of-date';
      else
         status.svn_oodmessage = '';
      end

      % next fields wc revision number (ended by name of file plus path ... might contain spaces, but is known)
      status.svn_rev = sscanf(reply{1}(9:end-length(file_fullpath)), '%s %*s %*s');
      
      % check if status.svn_status contains any invalid characters and if status.svn_rev is really a number
      validcharsfound = textscan(status.svn_status, validstatuschars, 'Whitespace', '');
      validcharsfound = validcharsfound{:};

      if isempty(validcharsfound) || ( cellfun(@length, validcharsfound) ~= length(status.svn_status) ) ||...
            (isnan(str2double(status.svn_rev)) && ~(status.svn_rev=='-'))
         for i = 1 : length(reply)
            disp(reply{i})
         end
         error('Unable to determine valid status/revision from output of SVN (output is listed above).');
      end

      % if svn found more than one matching file (should not be possible => just to make sure)
   elseif length(reply) > 2
      error('SVN found more than one file "%s"', status.filename);
   end
end

% add status to subversion cache if enabled
if globalsettings.core.svn_en_cache
   svn_cache.(fileext).(filename) = status;
end

