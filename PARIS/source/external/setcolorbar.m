% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates a colorbar ("EastOutside") in CurrentAxes and tags it as "colorbar"
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
% function setcolorbar()
%    Creates a colorbar in CurrentAxes and tags it as 'colorbar' to make it identifiable by savefigure.
% function setcolorbar(labels)
%    Creates a colorbar in CurrentAxes with user-defined YTickLabels given in LABELS (cell array of 
%    strings), and tags the colorbar as 'colorbar' to make it identifiable by savefigure.
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
% - colorbar with text on top (colorbar label)
%
% *******************************************************************************************************

function handle = setcolorbar(ticklabels)

% internal settings
internalsettings.fontsize = 10;
% internalsettings.fontname = 'Arial';
% internalsettings.fontname = 'Times New Roman';
internalsettings.fontname = 'Helvetica'; % Matlab default

% create colorbar
if nargin == 0
   handle = colorbar('EastOutside', 'FontSize',internalsettings.fontsize, 'FontName',internalsettings.fontname);
else
   handle = colorbar('EastOutside', 'FontSize',internalsettings.fontsize, 'FontName',internalsettings.fontname,...
      'YTickLabel',ticklabels);
end

% tag as colorbar
set(handle, 'Tag', 'colorbar');

