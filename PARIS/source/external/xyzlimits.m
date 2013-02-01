% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% calculate axes limits to fit data (plus small margins)
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
% lim = xyzlimits(data1, data2, data3, ...)
% lim = xyzlimits(data1, data2, data3, ..., 'linear')
%    Calculates axis limits to fit DATA1, DATA2, ... and returns the limits in LIM. 
%    Usage (e.g.) xlim(xyzlimit(data1, data2, data3, ...))
% lim = xyzlimits(data1, data2, data3, ..., 'log')
%    Calculates axis limits to fit DATA1, DATA2, ... and returns the limits for a log axis in LIM. 
%    Usage (e.g.) xlim(xyzlimit(data1, data2, data3, ..., 'log')) 
% 
%
%
% ***** Interface definition *****
% function lim = xyzlimits(varargin)
%    varargin   data vectors/matrices to fit axis to (x, y or z)
%               last entry (optional): 'log' for log axes, 'linear' for linear axes (default)
%    
%    lim        calculated axis limits
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
% - support for outliers (outlier filter optional, default off)
%
% *******************************************************************************************************

function lim = xyzlimits(varargin)

% internal settings
internalsettings.margin = 0.02;

% axis format (linear/log) given?
if ischar(varargin{end})
   nargin_max = nargin - 1;
   mode = varargin{end};
else
   nargin_max = nargin;
   mode = 'linear';
end

% reshape data to column vector
data = [];
for i = 1 : nargin_max
   if iscell(varargin{i})
      for j = 1 : numel(varargin{i})
         data = [data; varargin{i}{j}(:)]; %#ok<AGROW>
      end
   else
      data = [data; varargin{i}(:)]; %#ok<AGROW>
   end
end

% get minimum and maximum value of data
axes_max = max(data);
axes_min = min(data);

% check for max == min ~= 0
if axes_max <= axes_min
	lim = [axes_min ^    (1 + internalsettings.margin), ...
          axes_max ^ (1/(1 + internalsettings.margin))];
   lim = sort(lim); 
   return
end

% check for max == min == 0
if (axes_max == 0) && (axes_min == 0)
   lim = [-1 - internalsettings.margin,...
           1 + internalsettings.margin];
   return
end

% determine bounds for axis
switch(lower(mode))
   case 'linear'
      lim = [axes_min - (axes_max - axes_min) * internalsettings.margin,...
             axes_max + (axes_max - axes_min) * internalsettings.margin];
   case 'log'
      lim = [axes_min / (axes_max / axes_min) ^ internalsettings.margin,...
             axes_max * (axes_max / axes_min) ^ internalsettings.margin];
   otherwise
      warning('xyzlimits:UnknownMode',...
         'Unknown mode "%s" for axis, switching to linear. Valid options: linear, log.', lower(mode));
      lim = [axes_min - (axes_max - axes_min) * internalsettings.margin,...
             axes_max + (axes_max - axes_min) * internalsettings.margin];
end
