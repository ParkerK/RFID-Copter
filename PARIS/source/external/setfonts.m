% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% set all fonts in current figure (... Matlab does not save all fonts in a correct way)
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
% setfonts(fig)
%    Set all fonts in FIG to a default value ('latin modern sans')
% setfonts(fig, font)
%    Set all fonts in FIG to FONT.
% 
%
%
% ***** Interface definition *****
% function setfonts(fig, font)
%    fig    figure handle
%    font   (optional) font name; default: 'latin modern sans'
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

function setfonts(fig, font)

% input parameter checks / prepare input parameters
if nargin == 1
   font = 'latin modern sans'; % Matlab default would be Helvetica; that does not work on all systems
end

% set fonts
for h1 = findall(fig, 'Type', 'Axes')
   % axes text
   set(h1, 'FontName', font);
   % other text objects
   for h2 = findall(h1, 'Type', 'Text')
      set(h2, 'FontName', font);
   end
end
