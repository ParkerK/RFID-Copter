% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% logging system - version (output) system
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
% version_system() called from the command line
%    Outputs svn revision and content status (if svn is available) of all functions in
%    the current directory (pwd) and subdirectories to the command window.
% version_system() or version_system(print) called from a script
%    Outputs its own (= of version_system.m) version along with svn revision and content status  
%    (if svn is available) to the command window. Also adds the version of svn if it is available. PRINT
%    can be any string.
% version_system() called from a function
%    Outputs the version and (if svn is available) the svn version and content status of the calling 
%    function to the command window.
% version_system('all') called from a script
%    Like version_system() called from the command line
%
% Can be activated/deactivated via globalsettings.logging.versions (default: active)
%
%
% ***** Global Variables *****
% globalsettings
%    .logging.versions   activate version output (true/false)?
%    .logging.en_svn     include subversion revision information
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
% = extract version from header, not variable (will also fix some limitations of this fcn)
% ? redesign without dbstack() or get_stacklevel()
%
% *******************************************************************************************************

function version_system(print)
version = 'beta 3.0'; % version_system is the only core function that contains a version


% *******************************************************************************************************
% input parameter checks / prepare input parameters

if nargin == 0
   print = '';
end

   
% *******************************************************************************************************
% settings and information structs

% global settings
global globalsettings;

% check contents of globalsettings
%     prepare required data
expected.name = 'globalsettings';
expected.reqfields = {'logging', {'versions', 'en_svn'}};
%     check
errortext = contentcheck(globalsettings, expected);
%     output
if ~isempty(errortext)
   error('Incomplete globalsettings\n%s', errortext);
end

% check if version system should be active
if ~globalsettings.logging.versions   
   return
end

% necessary information about this function (svn fields will be completed later)
thisfcn = struct('name',mfilename, 'type', 'function', 'version',version, 'iscore', 1,...
                 'file',strcat(mfilename('fullpath'),'.m'), 'svn_isok',0, 'svn_version','', 'svn_status', '',...
                 'svn_rev','', 'svn_oodmessage','', 'svn_notcontrolled',0, 'notfound',0);

% necessary information about the calling function (will be completed later)
callingfcn = struct('level',0, 'name','', 'type', '', 'version','', 'iscore', 0, 'file','',...
                    'svn_isok',0, 'svn_version','', 'svn_status','', 'svn_rev','', 'svn_oodmessage','',...
                    'svn_notcontrolled',0,'notfound',0);


% *******************************************************************************************************
% gather information about the calling function

% get stack
stack = dbstack;
stack(1) = []; % remove this function (dbstack(1) does not really work with Matlab < 2007a)

% if this function was called from commandline: get version information of all files in pwd recursively
if isempty(stack)
   selfcalled(globalsettings.logging.en_svn, thisfcn, 'all');
   return
end

% we've been called from a function or script => get handle and function info
fhandle = str2func(stack(1).file(1:end-2));
finfo   = functions(fhandle);

% get stack level and functin name
callingfcn.level  = get_stacklevel(1);
callingfcn.name   = stack(1).file;
callingfcn.file   = finfo.file;

% get version and type information
callingfcn = mfile_status(callingfcn.name, callingfcn);

% if called from a script
if strcmpi(callingfcn.type, 'script')
   if strcmpi(print, 'all')
      selfcalled(globalsettings.logging.en_svn, thisfcn, 'all');
   else
      selfcalled(globalsettings.logging.en_svn, thisfcn, 'self');
   end
   return
end

% get information from svn (if available)
if globalsettings.logging.en_svn
   callingfcn = svn_status(callingfcn.file, callingfcn);
end


% *******************************************************************************************************
% output

output(globalsettings.logging.en_svn, callingfcn, '...');





% *******************************************************************************************************
% *******************************************************************************************************
% Operations if the primary function is called from script or command line
function selfcalled(en_svn, thisfcn, mode)

% get svn revision number and status for this function from svn (if enabled)
if en_svn; thisfcn = svn_status(thisfcn.file, thisfcn); end

switch lower(mode)
% own version and svn availability
   case 'self'
      % version
      if en_svn && thisfcn.svn_isok
         disp(sprintf('This is %s (rev.: %s | svn-rev.: %s | svn-status: "%s"%s)',...
            thisfcn.name, thisfcn.version, thisfcn.svn_rev, thisfcn.svn_status(1), thisfcn.svn_oodmessage));
      else
         disp(sprintf('This is %s (rev.: %s)', thisfcn.name, thisfcn.version));
      end
      % svn availability
      if thisfcn.svn_isok
         disp(sprintf('   ... SVN is available (version %s)', thisfcn.svn_version));
      else
         disp('   ... SVN is not available or not functional');
      end    
 
% version of all functions in pwd and subdirectories
   case 'all'
      % first of all: output our own version
      selfcalled(en_svn, thisfcn, 'self')
      % get a list of all matlab files recursively 
      files = rdir(fullfile(pwd, '**/*.m'));
      files = {files(:).name}'; % we're only interested in the names
      % get file names without path
      [paths, names] = cellfun(@fileparts, files, 'UniformOutput', false);
      % sort by names, index files accordingly
      [names, sortindices] = sort(names);
      % get maximum file name length (for spacer)
      maxlength = max(cellfun(@length, names));
      % thisfcn is used as template for loop: initialize fields that might not be set later
      fcn = thisfcn;
      fcn.file            = ''; % not used in loop
      fcn.level           = 0;    
      % get and print status of all files
      disp(sprintf('\nPrinting status of all .m-files in %s\nand subdirectories sorted by filename (%i files total):', pwd, length(names)));
      for i = 1 : length(names)
         % get necessary information from matlab        
         fcn.name = strcat(names{i}, '.m');
         fcn = mfile_status(files{sortindices(i)}, fcn);
         % get information from svn
         if en_svn; fcn = svn_status(files{sortindices(i)}, fcn); end
         % output
         if ~en_svn || ~thisfcn.svn_isok
            disp(sprintf('   %s%s(rev.: %18s)',...
               fcn.name, blanks(maxlength-length(fcn.name)+3), fcn.version));
%          elseif ~thisfcn.svn_isok
%             disp(sprintf('   %s%s(rev.: %18s |     svn is not functional/available     )',...
%                fcn.name, blanks(maxlength-length(fcn.name)+3), fcn.version));
         elseif fcn.notfound
            disp(sprintf('   %s%s(rev.: %18s |      SVN CANNOT FIND THIS FILE          )',...
               fcn.name, blanks(maxlength-length(fcn.name)+3), fcn.version));
         elseif fcn.svn_notcontrolled
            disp(sprintf('   %s%s(rev.: %18s |    file is not under version control    )',...
               fcn.name, blanks(maxlength-length(fcn.name)+3), fcn.version));
         else
            disp(sprintf('   %s%s(rev.: %18s | svn-rev.: %5s | svn-status: "%8s"%s)',...
               fcn.name, blanks(maxlength-length(fcn.name)+3), fcn.version,...
               fcn.svn_rev, fcn.svn_status, fcn.svn_oodmessage));
         end
      end
end





% *******************************************************************************************************
% *******************************************************************************************************
% status output (status: only first character = content status)
function output(en_svn, callingfcn, prefix)

% no SVN requested
if ~en_svn
   disp(sprintf('%s%s%s (rev.: %s)',...
      blanks(3*callingfcn.level), prefix, callingfcn.name, callingfcn.version));
% if svn is not available
elseif ~callingfcn.svn_isok
   disp(sprintf('%s%s%s (rev.: %s | svn not functional/available)',...
      blanks(3*callingfcn.level), prefix, callingfcn.name, callingfcn.version));
% if the file can't be found
elseif callingfcn.notfound
   disp(sprintf('%s%s%s (rev.: %s | SVN CANNOT FIND THIS FILE)',...
      blanks(3*callingfcn.level), prefix, callingfcn.name, callingfcn.version));
% if the file is not under version control
elseif callingfcn.svn_notcontrolled
   disp(sprintf('%s%s%s (rev.: %s | file is not under version control)',...
      blanks(3*callingfcn.level), prefix, callingfcn.name, callingfcn.version));
% if the file has not been modified and is not out-of-date
elseif callingfcn.svn_status(1) == ' ' && isempty(callingfcn.svn_oodmessage)
   disp(sprintf('%s%s%s (rev.: %s | svn-rev.: %s)', blanks(3*callingfcn.level),...
       prefix, callingfcn.name, callingfcn.version, callingfcn.svn_rev));
% if the file has been modified or is out-of-date
else
   disp(sprintf('%s%s%s (rev.: %s | svn-rev.: %s | svn-status: %s%s)', blanks(3*callingfcn.level),...
       prefix, callingfcn.name, callingfcn.version, callingfcn.svn_rev, callingfcn.svn_status(1),...
       callingfcn.svn_oodmessage));
end
