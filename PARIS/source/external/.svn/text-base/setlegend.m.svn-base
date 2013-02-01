% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates a legend in CurrentAxes and tags it as "legend"
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
% setlegend(legendtext)
%    like setlegend(legendtext, 'Best')
% setlegend(legendtext, location)
%    Creates a legend in CurrentAxes and tags is as 'legend' to make it identifiable by savefigure.
%    All strings in LEGENDTEXT are interpreted as TeX strings.  
% 
%
%
% ***** Interface definition *****
% function setlegend(legendtext, location)
%    legendtext   cell array of strings (>1 entry) or string (1 entry) for legend
%    location     (optional) position of legend; default: 'Best'
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

function setlegend(legendtext, location)

% input parameter checks / prepare input parameters
if nargin < 1
    error('Not enough input arguments');
elseif nargin == 1
    location = 'Best';
end

% internal settings
internalsettings.interpreter = 'tex';
internalsettings.fontsize = 10;
% internalsettings.fontname = 'Times New Roman';
internalsettings.fontname = 'Helvetica'; % Matlab default
% internalsettings.fontname = 'latin modern sans'; % looks better on Ubuntu

% create legend
handle = legend(legendtext,...
   'Location', location,...
   'Interpreter', internalsettings.interpreter,...
   'FontSize', internalsettings.fontsize,...
   'FontName', internalsettings.fontname);

% tag as legend
set(handle, 'Tag', 'legend');
