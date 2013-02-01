% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% call a matlab function to get its status
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
% status = mfile_status(file)
%    Tries to get the version and type (script, function) of file and returns the
%    information in a struct.
% status = mfile_status(file, status)
%    Tries to get the version and type (script, function) of file and adds the
%    information to struct STATUS (does only change its own fields of the struct).
%
%
% ***** Interface definition *****
% function status = mfile_status(file, status)
%    file          filename of file to check (with or without path and/or .m)
%    status        (optional) struct containing at least the fields stated below
%       .type      status returned by svn (first 7 characters)
%       .iscore    current working revision of the file
%       .version   a string like 'out of date' if file is out-of-date; empty string if not
%
%    status            struct containing the information obtained
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

function status = mfile_status(file, status)

% *******************************************************************************************************
% input parameters

% struct status
status_requiredfields = {'type', 'iscore', 'version'};
status_init   = {'', 0, ''};

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
% gather information about file

% split up path, filename, extension
[path, filename] = fileparts(file);

% get handle
fhandle = str2func(filename);

% find out if this file is a core function
status.iscore = iscorefunction(filename);

% get version information for non-corefunctions (there are no core-scripts) and find out if file is a
% function or a script
status.type = 'function';
if ~status.iscore
   try
      status.version = fhandle();
   catch %#ok<CTCH> % if the "function" is in fact not a function but a script
      status.type = 'script';
      status.version = '[script]';
   end
else
   status.version = '[core fcn]';
end
