% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates labels (x,y,z) and title in CurrentAxes
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
% setlabels(ttext, xtext, ytext, ztext)
%    Creates labels (x, y, optional z) and title in CurrentAxes. All strings are interpreted as TeX
%    strings. Use '' to skip the corresponding label/title.
% 
%
%
% ***** Interface definition *****
% function setlabels(ttext, xtext, ytext, ztext)
%    ttext   string for title
%    xtext   string for x label
%    ytext   string for y label
%    ztext   (optional) string for z label
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

function setlabels(ttext, xtext, ytext, ztext)

% internal settings
internalsettings.interpreter = 'tex';
internalsettings.titlesize = 11;
internalsettings.labelsize = 10;
internalsettings.fontname = 'Helvetica'; % Matlab default
% internalsettings.fontname = 'latin modern sans'; % looks better on Ubuntu

% input parameter checks / prepare input parameters
if nargin < 3
    error('Not enough input arguments');
elseif nargin == 3
    ztext = ' ';
end

% get CurrentAxes (create if necessary)
axeshandle = gca;

% set labels
if ~isempty(ttext)
   title(axeshandle, ttext,...
      'Interpreter', internalsettings.interpreter,...
      'FontSize', internalsettings.titlesize,...
      'FontName', internalsettings.fontname);
end
if ~isempty(xtext)
   xlabel(axeshandle, xtext,...
      'Interpreter', internalsettings.interpreter,...
      'FontSize', internalsettings.labelsize,...
      'FontName', internalsettings.fontname);
end
if ~isempty(ytext)
   ylabel(axeshandle, ytext,...
      'Interpreter', internalsettings.interpreter,...
      'FontSize', internalsettings.labelsize,...
      'FontName', internalsettings.fontname);
end
if ~isempty(ztext)
   zlabel(axeshandle, ztext,...
      'Interpreter', internalsettings.interpreter,...
      'FontSize', internalsettings.labelsize,...
      'FontName', internalsettings.fontname);
end



