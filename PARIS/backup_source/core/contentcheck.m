% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% check contents of a struct (e.g. check for complete settings)
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
% errortext = contentcheck(datastruct, expected)
%    Checks if DATASTRUCT contains all fields listed in EXPECTED.REQFIELDS. Returns an error message
%    listing all missing fields if at least one field is missing (additional fields are no problem). This
%    string can be used by the calling function to throw an error, making sure a stack trace shows a
%    problem in the calling function, not this function. Also supports structs in structs.
%    Note that this function consumes a lot of time for large cell arrays of structs; it can be switched
%    off by setting globalsettings.core.settings_check=false; 
%
%    Example with substruct:
%       settings = struct('main', 0, 'master', struct('sub',0)); % given
%       expected.name = 'settings struct';
%       expected.reqfields = {'main', 'master', {'sub'}}; % name of struct field before nested cellarray!
%       errortext = contentcheck(settings, expected);
%
%
% ***** Global Variables *****
% globalsettings
%    .core.settings_check    switch this function on?
%
%
% ***** Interface definition *****
% errortext = function contentcheck(datastruct, expected)
%    datastruct   struct to check (e.g. settings struct)
%    expected     verification info struct
%       .name        name of datastruct (e.g. 'settings'); displayed with error message
%       .reqfields   cell array with names of required fields (strings); use nested cell arrays for 
%                    substructs ({'fieldname of substruct', {'sub1', 'sub2}}); see example above
%
%    errortext   error text (empty if no errors were found); one line per error, grouped per substruct
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
% ? data type check needed
%
% *******************************************************************************************************

function errortext = contentcheck(datastruct, expected)

% *******************************************************************************************************
% input parameter checks / prepare paramters

% switched off?
global globalsettings;
if ~globalsettings.core.settings_check
   errortext = '';
   return
end

% is data a struct?
if ~isstruct(datastruct)
   errortext = sprintf('   Variable "%s" is not a struct or does not exist.', expected.name);
   return
end

% find cell arrays (=> substructs in datastruct)
subcells = cellfun(@iscell, expected.reqfields);
ind_subcells = find(subcells);

% check if there is always at least non noncell entry between cells
if diff(ind_subcells) <= 1
   criterr('expected.reqfields layout does not meet requirements. Check setup in calling function.');
end

% initialize errortext
errortext = '';


% *******************************************************************************************************
% check the first level (not a substruct)

% are all fields present?
if ~all(isfield(datastruct, expected.reqfields(~subcells)))
   % no, find out which ones are missing (more than expected is okay)
   names = expected.reqfields(~subcells);
   missing = names(~isfield(datastruct, expected.reqfields(~subcells)));
   % create error text; will lead to one line per missing field (%%s => %s, %s => name)
   %     template
   if isfield(expected, 'recursive') % this is a recursive call
      errortext = sprintf('   Required field "%%s" missing in substruct "%s".\n', expected.name);
   else
      errortext = sprintf('   Required field "%%s" missing in struct "%s".\n', expected.name);
   end
   %     errors
   errortext = sprintf(errortext, missing{:});
   %     remove last '\n'
   errortext = errortext(1:end-1);
end


% *******************************************************************************************************
% check substructs recursively; data ordering {main1, {sub1}, main2, {sub2, {subsub1}}, ...}

for i = 1 : sum(subcells)
   % prepare new expected struct
   subexp.name = expected.reqfields{ind_subcells(i) - 1};
   subexp.reqfields = expected.reqfields{ind_subcells(i)};
   subexp.recursive = true; % does not have to be boolean; the pure existence is enough
   % check if this field exists; skip recursion if it does not
   if ~isfield(datastruct, subexp.name)
      continue;
   end
   % recursive call
   subtxt = contentcheck(datastruct.(subexp.name), subexp);
   if ~isempty(subtxt)
      if isempty(errortext) % don't add unnecessary line breaks
         errortext = subtxt; 
      else
         errortext = sprintf('%s\n%s', errortext, subtxt);
      end
   end
end
