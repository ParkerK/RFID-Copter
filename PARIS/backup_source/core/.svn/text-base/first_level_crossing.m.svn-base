% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% determine first level crossing of a sampled function
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
% xl = first_level_crossing(x, y, yl)
%    Returns the first X where Y(X) == YL, NaN if no intersection was found. Tries to increase resolution 
%    using interpolation. 
%
%
% ***** Interface definition *****
% xl = first_level_crossing(x, y, yl)
%    x     x-values for function y(x); strictly monotonical
%    y     y-values for function y(x)
%    yl    level (y-value)
%    res   minimum resolution (quantization) of xl using spline interpolation
%
%    xl   intersection of y(x) and yl (x-value); NaN if no intersection was found
%   
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

function xl = first_level_crossing(x, y, yl, res)


% *******************************************************************************************************
% internal settings

internalsettings.maxres = 1e3; % maximum increase of resolution around intersection
internalsettings.nvic   =   2; % consider this number of samples around the intersection for interpolation (has to be enough to fit a spline)


% *******************************************************************************************************
% input parameter checks

% flatten x and y
x = x(:);
y = y(:);

% check dimensions
if numel(x) ~= numel(y)
   err('Number of elements in X and Y of Y(X) do not match.');
end
if numel(yl) ~= 1
   err('Level YL has to be scalar.');
end

% x has to be strictly monotonical for the interpolation to work
if ( sum(diff(x) < 0) ~= numel(x)-1 ) && ( sum(diff(x) > 0) ~= numel(x)-1 )
   err('X has to be strictly monotonical.');
end


% *******************************************************************************************************
% "the function"

% find first intersection
ind = find(y <= yl, 1, 'first'); % find intersection

% refine (if possible)
if ~isempty(ind) && ind > 1
   % check resolution; modify if necessary
   if abs(diff( x(max(1, min(length(x), ind-[1,0]))) )) / res > internalsettings.maxres
      res = abs(diff( x(max(1, min(length(x), ind-[1,0]))) )) / internalsettings.maxres;
      critwarn('Resolution exceeds limits. Truncating to %e', res);
   end
   % spline interpolation around the intersection (enough data to fit spline)
   ind_interp = max(1, ind-internalsettings.nvic) : min(length(x), ind+internalsettings.nvic-1);
   xi = x(ind-1) : res : x(ind);
   if length(xi) < 2; xi = x(ind-1:ind); end % res < diff
   yi = interp1(x(ind_interp), y(ind_interp), xi, 'spline', 'extrap');
   % find intersection in interpolated data (linear interp. should'nt be a problem in case of ambiguities
   xl  = interp1(yi, xi, yl, 'linear');
   % check interpolation (just to make sure)
   if xl < min(x(ind-[1,0])) || xl > max(x(ind-[1,0]))
      err('Interpolation has failed.');
   end
else
   xl = NaN;
end

% % % debug
% % close all; figure; hold on;
% % plot(x, y, 'b');
% % plot(x, ones(size(x))*yl, 'r');
% % plot(xi, yi, 'g');
% % plot(xl, yl, 'ro');
% % hold off; grid on;
