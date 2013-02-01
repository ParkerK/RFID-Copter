% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - estimate peak rms amplitude (via windowed power detection and rms=sqrt(power))
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
% version = peakampl()
%    Just returns the version number (string).
% ampl = peakampl(data, settings)
%    Returns the peak RMS (root mean square) of DATA according to SETTINGS. Estimation is done via 
%    windowed (!) power detection. For large vectors uniform/nonuniform downsampling can be done for
%    better performance (at the expense of estimation accuracy of course).
%
%
% ***** Interface definition *****
% function rms = peakrms(signal, settings)
%    signal     signal vector
%    settings   struct containing settings
%       .nwin      length of windows in samples
%       .ol        overlapping of windows in percent
%       .dsmode    downsampling mode {'uniform', 'rand', otherwise: no downsampling}
%                     uniform: signal = signal(1:settings.dsf:end)
%                     rand   : like uniform, but with random step size (on avgerage settings.dsf)
%       .dsf       downsampling factor
%
%    rms        estimated maximum RMS value
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
% - merge with other estimators ?
%
% *******************************************************************************************************

function rms = peakrms(signal, settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   rms = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'nwin', 'ol', 'dsmode', 'dsf'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% settings.dsf integer?
if settings.dsf ~= round(settings.dsf)
   warn('Downsampling factor (settings.dsf) has to be integer. Rounding.');
   settings.dsf = round(settings.dsf);
end


% *******************************************************************************************************
% peak detection

% if downsampling
switch lower(settings.dsmode)
   
   % uniform downsampling
   case 'uniform'
      warn('Downsampling by factor %g without anti-aliasing.', settings.dsf);
      signal = signal( round([1:settings.dsf:length(signal)]) );
      settings.nwin = round(settings.nwin / settings.dsf); % adapt window length...
      settings.ol = ceil(settings.nwin * settings.ol/100) / settings.nwin * 100; % ...and overlapping
   
   % random downsampling
   case 'rand'
      warn('Random downsampling by factor %i without anti-aliasing.', settings.dsf);
      signal = signal( round([1:settings.dsf:length(signal)]' + [rand(fix(length(signal)/settings.dsf), 1) * settings.dsf; 0]) );
      settings.nwin = round(settings.nwin / settings.dsf); % adapt window length...
      settings.ol = ceil(settings.nwin * settings.ol/100) / settings.nwin * 100; % ...and overlapping
      
   % no downsampling
   otherwise
end

% get avg power
%     setup est_power
settings.est_power.n  = settings.nwin;
settings.est_power.ol = settings.ol;
%     estimate amplitude via power
rms = sqrt(max(est_power(signal, settings.est_power)));
