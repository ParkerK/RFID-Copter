% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% check if SVN (subversion) is available + get its version number
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
% [available, version] = svn_isavailable()
%    if SVN is available:     available = 1 and version = version number (string)
%    if SVN is not available: available = 0 and version = ''
%
%
% ***** Global Variables *****
% globalsettings
%    .core.svn_isav_susp   check if svn is able to report the version of this function ("be suspicious")
%                          ... used to avoid problems with older SVN versions (available, but unable to
%                              deal with repositories created by newer versions)
%
%
% ***** Interface definition *****
% function [available, version] = svn_isavailable()
%    available   1 if svn is available, 0 if not
%    version     version of SVN if available, empty string if not
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
% ? deactivate SVN if not available for some time
%
% *******************************************************************************************************

function [available, version] = svn_isavailable()

% global settings
global globalsettings;

% try to call svn
[status, reply] = system('svn --version --quiet');

% check if svn is able to report the version of this function ("be suspicious")
if (status==0) && globalsettings.core.svn_isav_susp
   [status, msg] = system(sprintf('svn status --show-updates --verbose --quiet --non-interactive %s',...
      which(mfilename)));
   if status~=0
      reply = sprintf('SVN is available (%s), but reported an error: "%s"', strtrim(reply), strtrim(msg));
   end
end  

% get information
available = (status == 0);
version   = strtrim(reply); % remove end-of-line from string
