% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - calculate distances and angles between all given positions
%
% Notations: 
%    azimuth AZ is measured in the xy-plane, positive AZ is rotation towards y (counterclockwise)
%    elevation EL is measured from the positive z axis (0 is on z, pi/2 is in xy)
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
% version = get_distances()
%    Just returns the version number (string).
% dist_12 = get_distances(pos_1, pos_2)
%    Returns distances and angles between all objects with positions given POS_1 and all objects with 
%    position POS_2. POS_1 and POS_2 can be matrices [x1, y1, z1; x2, y2, z2; ...] or cell arrays 
%    {[x1,y1,z1], [x2,y2,z2], ...}. If POS_1 contains the position N objects and POS_2 contains M
%    positions, the returned matrices are of size N x M.
%    Note that the function is not limited to 3-dimensional positions, albeit angles are only calculated
%    in case of two (azimuth) and three (azimuth, elevation) dimensions.
%  
%
%
% ***** Function definition *****
% function [dist_12, az_12, el_12] = get_distances(pos_1, pos_2)
%    pos_1     position of N objects, e.g. 3-dim [x1, y1, z1; x2, y2, z2; ...; xN, yN, zN]
%              or cell array of vectors ; empty otherwise
%    pos_2     position of M objects, e.g. 3-dim [x1, y1, z1; x2, y2, z2; ...; xM, yM, zM]
%              or cell array of vectors
%
%    dist_12   NxM distance matrix; distance between objects (rows: pos_1, columns: pos_2)
%    az_12     NxM matrix with azimuth angles [rad] between objects (only for 2,3-dim; empty otherwise)
%    el_12     NxM matrix with elevation angles [rad] between objects (only for 3-dim; empty otherwise)
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

function [dist_12, az_12, el_12] = get_distances(pos_1, pos_2)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   dist_12 = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameters

if nargin ~= 2
   criterr('Wrong number of input arguments.');
end

% empty matrices
if isempty(pos_1)
   warn('Call parameter pos_1 is empty => distance matrix will be empty.');
end
if isempty(pos_2)
   warn('Call parameter pos_2 is empty => distance matrix will be empty.');
end

% convert cell arrays to matrices if necessary
if iscell(pos_1)
   % make sure rows are readers ...
   pos_1 = pos_1(:);
   % ... and coordinates are column vectors ...
   pos_1 = cellfun(@(x) x(:)', pos_1, 'UniformOutput',false);
   % ... before transformation
   pos_1 = cell2mat(pos_1);
end
if iscell(pos_2)
   % make sure rows are readers ...
   pos_2 = pos_2(:);
   % ... and coordinates are column vectors ...
   pos_2 = cellfun(@(x) x(:)', pos_2, 'UniformOutput',false);
   % ... before transformation
   pos_2 = cell2mat(pos_2);
end

% make sure pos_1 and pos_2 have identical dimensionality
if size(pos_1, 2) ~= size(pos_2, 2)
   err('Dimensionality of pos_1 (%i) and pos_2 (%i) does not match.', size(pos_1, 2), size(pos_2, 2));
end

% *******************************************************************************************************
% calculate distances
% ... this function is not critical for performance => not vectorized for simplicity

% pre-allocate matrices
dist_12 = zeros(size(pos_1,1), size(pos_2,1));
switch size(pos_1, 2)
   case 2
      az_12 = zeros(size(pos_1,1), size(pos_2,1));
      el_12 = [];
   case 3
      az_12 = zeros(size(pos_1,1), size(pos_2,1));
      el_12 = zeros(size(pos_1,1), size(pos_2,1));
   otherwise
      az_12 = [];
      el_12 = [];
end

% calculate distances and angles
for i = 1 : size(pos_1, 1)
   for j = 1 : size(pos_2, 1)
      d_12 = pos_2(j,:) - pos_1(i,:); % vector 1->2
      switch length(d_12)
         case 3 % {distance, azimuth, elevation}
            [az_12(i,j), el_12(i,j), dist_12(i,j)] =...
               cart2sph(d_12(1), d_12(2), d_12(3));
            el_12(i,j) = pi/2 - el_12(i,j); % change notation
         case 2 % {distance, azimuth}
            [az_12(i,j), dist_12(i,j)] =...
               cart2pol(d_12(1), d_12(2));
         otherwise % {distance}
            dist_12(i,j) = norm(d_12);
      end
   end
end

