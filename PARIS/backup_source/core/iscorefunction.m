% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% determine if a function is part of the core functions (only primary functions)
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
% iscorefunction(function_handle)
% iscorefunction(file name plus path and/or extension)
% iscorefunction(function name)
%    returns 1 if the function is part of the core files
%
% ATTENTION: If no core path is provided in globalsettings.path.core, the function assumes that 
%            it resides within the core-root-directory (corepath = path of this function).
%            Core files are per definition files that reside in corepath and subdirectories.
%
%
% ***** Global Variables *****
% globalsettings
%    .path.corepath   root path for core files
%
%
% ***** Interface definition *****
% function itis = iscorefunction(fcn)
%    fcn    function handle / function name / file name of function to check
%
%    itis   1 if fcn is a core function, 0 if not
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

function itis = iscorefunction(fcn)

% *******************************************************************************************************
% settings

% global settings
global globalsettings;

% check if a core path has been provided ... if not: use pwd as corepath
if isstruct(globalsettings) && isstruct(globalsettings.path) && isfield(globalsettings.path, 'core')
   corepath = globalsettings.path.core;
else
   corepath = fileparts(mfilename('fullpath'));
end


% *******************************************************************************************************
% function

% if fcn is a function handle: get filename (without extension)
if isa(fcn, 'function_handle')
   finfo = functions(fcn);
   fcn = finfo.file;
end

% remove path + extension (they are unnecessary)
[path, filename] = fileparts(fcn);

% get a list of all files in corepath and subdirectories
corefiles = rdir(fullfile(corepath, '**/*.m'));
corefiles = {corefiles(:).name}'; % we're only interested in the names
[paths, corefiles] = cellfun(@fileparts, corefiles, 'UniformOutput', false);

% find out if filename is part of corefiles
itis = 0;
for i = 1 : length(corefiles)
   if strcmp(filename, corefiles{i}) == 1
      itis = 1;
      break;
   end
end
