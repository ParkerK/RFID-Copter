% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% initialize paths and global variables for ranging m-files
%
%        Note that there is hardcoded information in here!
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
% globalinit()
%    Returns without any action taken ... prevents execution by function VERSION_SYSTEM.
% globalinit('sim');
%    Defines paths and global variables and switches on all warnings and emails (no messages). SVN and
%    matfile caches are activated in this mode. Settings structs are not checked (takes too much time).
% globalinit('sim-silent');
%    Defines paths and global variables and switches off all messages, warnings, and emails. All caches
%    (SVN, matfile-loading) are activated in this mode. Settings structs are not checked (takes too much
%    time).
% globalinit('sim-verbose');
%    Defines paths and global variables and switches on all messages, warnings, and emails. SVN and
%    matfile caches are activated in this mode. Settings structs are not checked (takes too much time).
% globalinit('silent');
%    Defines paths and global variables and switches off all messages, warnings, and emails. All caches
%    (SVN, matfile-loading) are deactivated. Settings structs are checked.
% globalinit('exceptions');
%    Defines paths and global variables. Only error and warning messages will be displayed. Emails are
%    active. All caches (SVN, matfile-loading) are deactivated. Settings structs are checked.
% globalinit('verbose');
%    Defines paths and global variables and switches on all messages, warnings and emails. All caches
%    (SVN, matfile-loading) are deactivated. Settings structs are checked.
% globalinit('UC');
%    Defines paths and global variables and switches on all messages and warnings. No emails will be
%    sent, the simulator will enter debug mode in case of an error. All caches (SVN, matfile-loading)
%    are deactivated. Settings structs are checked.
% globalinit(string);
%    Defines paths and global variables and switches on all warnings; messages will not be displayed. 
%    The email system is on. STRING can be any string (e.g., 'init'). All caches (SVN, matfile-loading)
%    are deactivated. Settings structs are checked.
%
%
% ***** Global Variables *****
% globalsettings
%    .path                  path strings
%       .master                root path (:= path to this function)
%       .core                  path to core functions
%       .simulator             path to simulator root
%    .logging               settings related to version/message (output) system
%       .versions              enable/disable version outputs
%       .en_svn                (true) add/do not add subversion revision info in version output
%       .exceptions            enable/disable error / critical error outputs 
%                              enable/disable stopping of diary in case of an error
%       .warnings              enable/disable warning / critical warning outputs
%       .headlines             enable/disable headline outputs in main modules
%       .messages              enable/disable message outputs
%    .core                  settings for core functions
%       .svn_isav_susp         (false) suspicious availability detection for SVN; also checks for version
%                              incompatibilities 
%       .svn_head              (false) include svn server in svn calls (compare to HEAD; server load!)
%       .svn_en_cache          (true) call subversion only once per file; store result in cache 
%       .mat_en_cache          (true) load matfiles only once; store in cache for future loads
%       .mat_stubborn          ignore missing files; return empty struct if file is missing
%       .settings_check        check settings (contentcheck.m)? ...takes a lot of time for large settings
%       .mailto                email address for simulator messages
%       .subj_prefix           subject prefix for email identification (allows for automatic filtering)
%       .debug                 start in debug mode (special options within some functions)
%       .sfir_char             filename for sparse FIR performance characteristic
%       .ticket_timeframe      only one ticket-renewal per timeframe (minutes)
%       .ticket_retries        number of immediate retries for failed renewals
%       .ticket_waittime       seconds waittime between immediate retries for failed renewals (seconds)
%       .ticket_tolerance      failed renewals within this tolerance time will not be reported (minutes)
%       .ticket_minrenewleft   minimum renewal time left before user interaction is requested (minutes)
%    
%    .misc              miscellaneous settings / information
%       .pid                process ID of Matlab (if Unix's screen is used: PID of screen is .pid-1)
%       .hostname           name of host we're currently running on
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

function settings = globalinit(command)

% prevent execution by VERSION_SYSTEM
if nargin == 0 || ~ischar(command)
   return
end

% global struct
clear global globalsettings;
global globalsettings;


% *******************************************************************************************************
% master path

% "master" path = path to this file
globalsettings.path.master = fileparts(mfilename('fullpath'));
% find all subdirectories of "master" path that are not .svn directories
%     split up (this might take a few seconds)
allpaths = textscan(genpath(globalsettings.path.master), '%s', 'delimiter', pathsep);
allpaths = allpaths{:};
%     collect everything that does not contain '.svn', resolve symbolic links
goodpaths = '';
for i = 1 : length(allpaths)
   if isempty(findstr(allpaths{i}, '.svn'))
      % let the system resolve symbolic links on Unix systems (otherwise Debugger might not work properly)
      if isunix
         [status, newpath] = system(sprintf('readlink -f "%s"', allpaths{i}));
      else
         newpath = allpaths{i};
      end
      
      % add to path list suited for addpath
      goodpaths = strcat(goodpaths, newpath, pathsep);
   end
end
% add to Matlab path (reset old path first)
restoredefaultpath;
addpath(goodpaths);


% *******************************************************************************************************
% other paths

% core files (directory marked by a file)
if isunix
   [status, reply] = system(sprintf('find %s -name ".SPR_COREPATH"', globalsettings.path.master));
   markerfile = textscan(reply, '%s', 'delimiter', '\n', 'Whitespace', '');
   markerfile = markerfile{1};
else
   markerfile = rdir(fullfile(globalsettings.path.master, '**/.SPR_COREPATH')); % quite slow
   markerfile = {markerfile(:).name}'; % we're only interested in the name
end
if length(markerfile) == 1
   globalsettings.path.core = fileparts(markerfile{1});
else
   error('Could not find Core path or found more than one. Check markerfile ".SPR_COREPATH".');
end

% simulator files (base directory marked by a file)
if isunix
   [status, reply] = system(sprintf('find %s -name ".SPR_SIMPATH"', globalsettings.path.master));
   markerfile = textscan(reply, '%s', 'delimiter', '\n', 'Whitespace', '');
   markerfile = markerfile{1};
else
   markerfile = rdir(fullfile(globalsettings.path.master, '**/.SPR_SIMPATH')); % quite slow
   markerfile = {markerfile(:).name}'; % we're only interested in the name
end
if length(markerfile) == 1
   globalsettings.path.simulator = fileparts(markerfile{1});
else
   error('Could not find Simulator path or found more than one. Check markerfile ".SPR_SIMPATH".');
end


% *******************************************************************************************************
% logging system
switch lower(command)
   % do not print anything
   case {'silent', 'sim-silent'}
      globalsettings.logging.versions   = 0;
      globalsettings.logging.en_svn     = 0;
      globalsettings.logging.exceptions = 0;
      globalsettings.logging.warnings   = 0;
      globalsettings.logging.headlines  = 0;
      globalsettings.logging.messages   = 0;
      
   % print only exceptions (warnings, errors)
   case {'exceptions', 'sim'}
      globalsettings.logging.versions   = 0;
      globalsettings.logging.en_svn     = 0;
      globalsettings.logging.exceptions = 1;
      globalsettings.logging.warnings   = 1;
      globalsettings.logging.headlines  = 1;
      globalsettings.logging.messages   = 0;
      
   % print everything
   case {'verbose', 'uc', 'sim-verbose'}
      globalsettings.logging.versions   = 1;
      globalsettings.logging.en_svn     = 1;
      globalsettings.logging.exceptions = 1;
      globalsettings.logging.warnings   = 1;
      globalsettings.logging.headlines  = 1;
      globalsettings.logging.messages   = 1;

   % standard setting: activate everything except messages
   otherwise
      globalsettings.logging.versions   = 1;
      globalsettings.logging.en_svn     = 1;
      globalsettings.logging.exceptions = 1;
      globalsettings.logging.warnings   = 1;
      globalsettings.logging.headlines  = 1;
      globalsettings.logging.messages   = 0;
end


% *******************************************************************************************************
% exception handling, email system
switch lower(command)
   % simulation mode: send error messages per mail, enter debug mode in case of an error
   case {'sim', 'sim-verbose'}
      globalsettings.core.mailto      = '';
      globalsettings.core.subj_prefix = '[SIM] '; % prefix for subject to identify email in clients
      globalsettings.core.debug       = false;
         
   % "under construction" mode: do not send email messages, enter debug mode in case of an error
   case 'uc'
      globalsettings.core.mailto      =   '';
      globalsettings.core.subj_prefix =   '';
      globalsettings.core.debug       = true;
      
   % standard setting: deactivate exception handling and emails
   otherwise
      globalsettings.core.mailto      =    '';
      globalsettings.core.subj_prefix =    '';
      globalsettings.core.debug       = false;
end


% *******************************************************************************************************
% misc

% reset alls caches
clear global svn_cache; % subversion status
clear global loadmat_cache; % matfile loading
clear global ticket_cache; % ticket renewal

% determine Matlab's PID and the hostname we're running on
globalsettings.misc.pid = feature('getpid');
if isunix
   [dummy, globalsettings.misc.hostname] = system('echo $HOSTNAME');
   globalsettings.misc.hostname = deblank(globalsettings.misc.hostname);
else
   globalsettings.misc.hostname = '[unknown]';
end

% Matlab sets identical initial random states at each startup => randomize here
RandStream.setDefaultStream(RandStream('mt19937ar','seed',sum(100*clock)));

% set default graphics position in twoscreen mode (if not already done)
resolution = get(0, 'ScreenSize');
position = get(0, 'DefaultFigurePosition');
if resolution(3)/resolution(4) > 2 && position(1) + position(3) > resolution(3)/2
   position(1) = resolution(3)/4 - position(3)/2;
   set(0, 'DefaultFigurePosition', position);
end


% *******************************************************************************************************
% core functions

switch lower(command)
   % simulation modes: activate all caches
   case {'sim', 'sim-silent', 'sim-verbose'}
      % subversion
      globalsettings.core.svn_isav_susp = false; % be suspicious when checking availability of SVN
      globalsettings.core.svn_head      = false; % svn-status compared to HEAD (true) or BASE (false)
      globalsettings.core.svn_en_cache  =  true; % cache for subversion (only one call per file)
      
      % matfile-loading
      globalsettings.core.mat_en_cache =  true; % cache for matfile loading
      globalsettings.core.mat_stubborn = false; % "stubborn" matfile loading; load emtpy struct if file is missing
      
      % settings content check (takes time for large settings structs)
      globalsettings.core.settings_check = false; % do not check settings; simulator should work! (takes too long)
      
      % standard setting: deactivate all caches (might cause undesired behavior, e.g, modified files not loaded)
   otherwise
      % subversion
      globalsettings.core.svn_isav_susp = false; % be suspicious when checking availability of SVN
      globalsettings.core.svn_head      = false; % svn-status compared to HEAD (true) or BASE (false)
      globalsettings.core.svn_en_cache  = false; % cache for subversion (only one call per file)
      
      % matfile-loading
      globalsettings.core.mat_en_cache = false; % cache for matfile loading
      globalsettings.core.mat_stubborn = false; % "stubborn" matfile loading; load emtpy struct if file is missing
      
      % settings content check (takes time for large settings structs)
      globalsettings.core.settings_check = true; % check settings
end

% performance characteristics for sparse FIR implementation
globalsettings.core.sfir_char = sprintf('syschar_sparsefir_%s', globalsettings.misc.hostname);

% ticket renewal
%     only one renewal per timeframe
globalsettings.core.ticket_timeframe    =    60; % minutes
%     immediate retries for failed renewals
globalsettings.core.ticket_retries      =     1; % number of retries
globalsettings.core.ticket_waittime     =    30; % seconds waittime between renewals
%     failed renewals within this tolerance time will not be reported
globalsettings.core.ticket_tolerance    = 12*60; % minutes
%      minimum renewal time left before user interaction is requested
globalsettings.core.ticket_minrenewleft = 48*60; % minutes

% return a copy of these settings
settings = globalsettings;
   
