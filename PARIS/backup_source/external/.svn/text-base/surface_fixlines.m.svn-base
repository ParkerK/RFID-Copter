% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% workaround for bug in Matlab surface plots (rendering)
%    When viewing a surface with linear axes from above (i.e., view([0,90])), intersections with a flat
%    minimum (e.g. a saturated minimum) are visible as white lines (unwanted transparency at intersection 
%    of polygons).
%    Possible reference: http://www.opengl.org/resources/faq/technical/polygonoffset.htm
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
% surface_fixlines(x, y, z)
%    Insert a plane below min(min(Z)) spanning the entire range of X and Y.
% surface_fixlines(x, y, z)
%    Insert a plane over the entire range X, Y, Z; color min(c). Only for simple surfaces!
%
%
% ***** Interface definition *****
% function surface_fixlines(x, y, z)
%    x   vector/matrix containing x axis data
%    y   vector/matrix containing y axis data
%    z   vector/matrix containing z axis data / matrix spanning the surface z(x,y)
%    c   (optional)  vector/matrix containing color axis data for surface c(x,y,z)
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

function surface_fixlines(x, y, z, c)

if nargin == 3
   surface([min(min(x)), max(max(x))], [min(min(y)), max(max(y))], ones(2,2)*min(min(z))-eps, 'EdgeColor','none');
elseif nargin == 4
   x1 = min(x(:)); x2 = max(x(:));
   y1 = min(y(:)); y2 = max(y(:));
   z1 = min(z(:)); z2 = max(z(:));
   if abs(x2 - x1) < 1000*eps
      [Y,Z] = meshgrid([y1,y2], [z1,z2]);
      X = x1 * ones(2,2);
   elseif abs(y2 - y1) < 1000*eps
      [X,Z] = meshgrid([x1,x2], [z1,z2]);
      Y = y1 * ones(2,2);
   elseif abs(z2 - z1) < 1000*eps
      [X,Y] = meshgrid([x1,x2], [y1,y2]);
      Z = z1 * ones(2,2);
   else
      error('Won''t work for this surface...');
   end
   surface(X, Y, Z, ones(2,2)*min(c(:)), 'EdgeColor','none');
else
   error('Wrong number of input arguments.');
end
