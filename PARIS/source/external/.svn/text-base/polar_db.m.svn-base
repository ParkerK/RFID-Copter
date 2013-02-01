% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% polar plot in dB (very simple implementation still using Matlab's POLAR)
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
% polar_db(angle, gain, style)
%    Like POLAR(angle, gain, style), but supports negative values for gain (typ. gain patterns). The 
%    function relies on an accurate prediction of the labels Matlab's POLAR function produces, 
%    so it may fail for future implementations of POLAR (hopefully mathworks will be then have replaced
%    POLAR with something that can handle negative values). The dynamic range of the plot is fixed to
%    25 dB for simplicity.
%
%
% ***** Interface definition *****
% function polar_db(angle, gain, style)
%    angle   angle of polar plot in radians
%    gain    gain to plot in dB
%    style   linestyle
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

function polar_db(angle, gain, style)

% internal settings
internalsettings.round =  5; % dB rounding interval of dB labels
internalsettings.scale = 25; % dB dynamic range

% calculate dB scale
dblabels = round( max(gain) / internalsettings.round ) * internalsettings.round + ...
   - internalsettings.scale : internalsettings.round : 0;

% scale and shift gain
gain = gain - max(gain) + internalsettings.scale;
gain(gain < 0) = 0;

% polar plot
polar(angle, gain, style);

% replace strings in dB scale
%     make sure the handles can't be hidden
old_shh = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');
%     replace
for i = 0 : round(internalsettings.scale / internalsettings.round)
   oldstr = ['  ', num2str(i * internalsettings.round)]; % two whitespaces to separate from angles
   newstr = [' ', num2str(dblabels(i+1)), ' dB']; % make that one ... we don't want to replace iteratively
   if ~isempty(findobj(gca, 'String', oldstr))
      set(findobj(gca, 'String', oldstr),'String', newstr);
   end
end
%     reset ShowHiddenHandles
set(0, 'ShowHiddenHandles', old_shh);
