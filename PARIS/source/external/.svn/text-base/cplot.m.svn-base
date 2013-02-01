% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% extended "custom" plot() command allowing easy and fast control of parameters
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
% Behaves similar to PLOT, but redefines certain settings (optional). Please note that no input 
% parameter checks are done. If only LINE is defined, MARKER will be set to 'none' and vice versa.
% If X and/or Y are matrices and LINE and/or MARKER are defined, all plots will be set according 
% to LINE and/or MARKER (and thus look the same).
%
% Usage:
%    keep everything standard (e.g. like plot(x)):
%       - cplot(y)
%       - cplot(x,y)
%    redefine only line settings:
%       - cplot(y, line)
%       - cplot(x, y, line}
%    redefine only marker settings:
%       - cplot(y, {}, marker}
%       - cplot(x, y, {}, marker)
%    redefine line and marker settings:
%       - cplot(y, line, marker)
%       - cplot(x, y, line, marker)
% 
%
%
% ***** Interface definition *****
% function cplot(x, y, line, marker)
%    x        (optional) x data for plot (scalar, vector, matrix)
%    y        y data for plot (scalar, vector, matrix)
%    line     (optional) cell array {LineStyle, LineWidth, Color}
%    marker   (optional) cell array {Marker, Markersize, MarkerEdgeColor=MarkerFaceColor}
%                                or {Marker, Markersize, MarkerEdgeColor, MarkerFaceColor}
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
% - input parameter checks ?
%
% *******************************************************************************************************

function cplot(varargin)


if nargin < 1
   error('Not enough input arguments');

   % plot(y)
elseif nargin == 1
   plot(varargin{1});

   % plot(x,y) or plot(y, {line})
elseif nargin == 2
   if iscell(varargin{2})
      plot(varargin{1}, ...
         'LineStyle',varargin{2}{1}, 'LineWidth', varargin{2}{2}, 'Color',varargin{2}{3}, ...
         'Marker', 'none');
   else
      plot(varargin{1}, varargin{2});
   end

   % plot(y, {line}, {marker}) or plot(y, {}, {marker}) or plot(x,y,{line})
elseif nargin == 3
   if iscell(varargin{2})
      if isempty(varargin{2})
         if length(varargin{3}) == 3
            plot(varargin{1},...
               'LineStyle','none', ...
               'Marker',varargin{3}{1}, 'Markersize', varargin{3}{2}, 'MarkerEdgeColor',varargin{3}{3}, 'MarkerFaceColor',varargin{3}{3});
         else
            plot(varargin{1},...
               'LineStyle','none', ...
               'Marker',varargin{3}{1}, 'Markersize', varargin{3}{2}, 'MarkerEdgeColor',varargin{3}{3}, 'MarkerFaceColor',varargin{3}{4});
         end
      else
         if length(varargin{3}) == 3
            plot(varargin{1},...
               'LineStyle',varargin{2}{1}, 'LineWidth', varargin{2}{2}, 'Color',varargin{2}{3}, ...
               'Marker',varargin{3}{1}, 'Markersize', varargin{3}{2}, 'MarkerEdgeColor',varargin{3}{3}, 'MarkerFaceColor',varargin{3}{3});
         else
            plot(varargin{1},...
               'LineStyle',varargin{2}{1}, 'LineWidth', varargin{2}{2}, 'Color',varargin{2}{3}, ...
               'Marker',varargin{3}{1}, 'Markersize', varargin{3}{2}, 'MarkerEdgeColor',varargin{3}{3}, 'MarkerFaceColor',varargin{3}{4});
         end
      end
   else
      plot(varargin{1}, varargin{2}, ...
         'LineStyle',varargin{3}{1}, 'LineWidth', varargin{3}{2}, 'Color',varargin{3}{3}, ...
         'Marker', 'none');
   end

   % plot(x,y,{},{marker}) or plot(x,y,{line},{marker})
elseif nargin == 4
   if isempty(varargin{3})
      if length(varargin{4}) == 3
      plot(varargin{1}, varargin{2}, ...
         'LineStyle', 'none',...
         'Marker',varargin{4}{1}, 'Markersize', varargin{4}{2}, 'MarkerEdgeColor',varargin{4}{3}, 'MarkerFaceColor',varargin{4}{3});
      else
         plot(varargin{1}, varargin{2}, ...
         'LineStyle', 'none',...
         'Marker',varargin{4}{1}, 'Markersize', varargin{4}{2}, 'MarkerEdgeColor',varargin{4}{3}, 'MarkerFaceColor',varargin{4}{4});
      end
   else
      if length(varargin{4}) == 3
         plot(varargin{1}, varargin{2}, ...
            'LineStyle',varargin{3}{1}, 'LineWidth', varargin{3}{2}, 'Color',varargin{3}{3},...
            'Marker',varargin{4}{1}, 'Markersize', varargin{4}{2}, 'MarkerEdgeColor',varargin{4}{3}, 'MarkerFaceColor',varargin{4}{3});
      else
         plot(varargin{1}, varargin{2}, ...
            'LineStyle',varargin{3}{1}, 'LineWidth', varargin{3}{2}, 'Color',varargin{3}{3},...
            'Marker',varargin{4}{1}, 'Markersize', varargin{4}{2}, 'MarkerEdgeColor',varargin{4}{3}, 'MarkerFaceColor',varargin{4}{4});
      end
   end

else
   error('Too many input arguments');
end
   
