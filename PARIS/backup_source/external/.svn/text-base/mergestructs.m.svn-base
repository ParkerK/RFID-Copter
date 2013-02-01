% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% merge two structs (e.g., settings)
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
% struct1 = mergestructs(struct1, struct2)
%    Merge STRUCT2 into STRUCT1, overwriting existing fields in STRUCT1.
% struct1 = mergestructs(struct1, struct2, fields)
%    Merge fields of STRUCT2 specified by FIELDS into STRUCT1, overwriting existing fields in STRUCT1.
%
%
% ***** Interface definition *****
% function struct1 = mergestructs(struct1, struct2, fields)
%    struct1   base struct
%    struct2   struct that should be merged with STRUCT1
%    fields    (optional) cell array of field names that should be merged
%
%    struct1   STRUCT1 merged with STRUCT2 
%
%
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% alpha 1.0  2012-05-17   arnitz      
%
%
% ***** Todo *****
%
% *******************************************************************************************************

function struct1 = mergestructs(struct1, struct2, fields)

% safety check: make sure both settings and optsettings are structs
if ~isstruct(struct1) || ~isstruct(struct2)
   error('Both STRUCT1 and STRUCT2 have to be structs.');
end

% names not specified => copy all settings
if nargin == 2
   fields = fieldnames(struct2);
end

for i = 1 : length(fields)
   % skip non-existing field names ...
   if ~isfield(struct2, fields{i})
      continue;
   else % ... otherwise: copy
      struct1.(fields{i}) = struct2.(fields{i});
   end
end

% sort field names
struct1 = orderfields(struct1);
