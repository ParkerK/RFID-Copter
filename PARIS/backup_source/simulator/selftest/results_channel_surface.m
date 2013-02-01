% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_channel_reflection function
% !!! EXISTING RESULTS WILL BE OVERWRITTEN - USE WITH CARE !!!
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


% *******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% initialize global stuff, print only warnings and errors (assumes that globalinit.m is part of path)
globalinit('exceptions');

% "switch to here"
cd(fileparts(mfilename('fullpath')));


% *******************************************************************************************************
% test setup

% output directory and filename
directory = 'results';
filename  = 'results_channel_surface';

% modes for surfaces (have to be identical to settings in CHANNEL_MAIN and CHANNEL_NEWPOS)
internalsettings.surf_off  = 0; % neither reflection, nor transmission => surface inactive
internalsettings.surf_tran = 1; % "reflect" mode
internalsettings.surf_refl = 2; % "transmit" mode
 
% random test setup bounds [min, max] (uniformly distributed)
testsettings.naoi   = [1, 10]; % number of gains to be calculated in parallel (size of .aoi NxM)
testsettings.n2     = [1, 10]; % refractive index of second medium (first: air)

% partitioning
partitioning.names = {'off, perp. pol.', 'off, par. pol.',...
   'refl., perp. pol.', 'refl., par. pol.', 'transm., perp. pol.', 'transm, par. pol.'};
partitioning.runs  = [100, 100, 1000, 1000, 1000, 1000];


% *******************************************************************************************************
% complete/check partitioning

if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check partitioning.');
end

% indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% create settings

for j = 1 : length(partitioning.names)
   for i = 1 + partitioning.indices(j) : partitioning.indices(j+1)
      settings{i}.aoi = 90 * rand(randi(testsettings.naoi), randi(testsettings.naoi)); %#ok<*SAGROW>
      settings{i}.n2  = testsettings.n2(1) + diff(testsettings.n2) * rand;
      settings{i}.dim = randi([1,3]); % "3D"
      % mode and polarization-dependent settings
      switch lower(partitioning.names{j})
         case 'off, perp. pol.' % off, perpendicular polarization
            results{i}.mode     = 'off'; % store independently to include internalsettings in check (may have changed)
            settings{i}.mode    = internalsettings.surf_off;
            settings{i}.pol_dim = settings{i}.dim; % surface normal dimension == polarization dimension
         case 'off, par. pol.' % off, parallel polarization
            results{i}.mode     = 'off';
            settings{i}.mode    = internalsettings.surf_off;
            settings{i}.pol_dim = 1 + mod(settings{i}.dim - 2, 2); % surface normal dimension ~= polarization dimension
         case 'refl., perp. pol.' % reflection, perpendicular polarization
            results{i}.mode     = 'reflection';
            settings{i}.mode    = internalsettings.surf_refl;
            settings{i}.pol_dim = settings{i}.dim;
         case 'refl., par. pol.' % reflection, parallel polarization
            results{i}.mode     = 'reflection';
            settings{i}.mode    = internalsettings.surf_refl;
            settings{i}.pol_dim = 1 + mod(settings{i}.dim - 2, 2);
         case 'transm., perp. pol.' % transmission, perpendicular polarization
            results{i}.mode     = 'transmission';
            settings{i}.mode    = internalsettings.surf_tran;
            settings{i}.pol_dim = settings{i}.dim;
         case 'transm, par. pol.' % transmission, parallel polarization
            results{i}.mode     = 'transmission';
            settings{i}.mode    = internalsettings.surf_tran;
            settings{i}.pol_dim = 1 + mod(settings{i}.dim - 2, 2);
         otherwise
            error('Unsupported partitioning.names "%s".', lower(partitioning.names{j}));
      end
   end
end


% *******************************************************************************************************
% calculate expected results (equations checked manually by creating a plot for aoi 0..90)

for j = 1 : length(partitioning.names)
   for i = 1 + partitioning.indices(j) : partitioning.indices(j+1)
      
      % extract n1, n2, and aoi for convenience
      n1  = 1;
      n2  = settings{i}.n2;
      aoi = settings{i}.aoi * pi/180;
      
      % check polarization / calculate gain factor (Fresnel equations)
      switch lower(partitioning.names{j}(end-9:end))
         case 'perp. pol.'
            if settings{i}.pol_dim ~= settings{i}.dim
               error('Polarization should be perpendicular, but normal vector directions does not match.');
            end
            results{i}.gain = (n1*sqrt(1-(n1/n2*sin(aoi)).^2)-n2*cos(aoi)) ./ (n1*sqrt(1-(n1/n2*sin(aoi)).^2)+n2*cos(aoi));
         case ' par. pol.'
            if settings{i}.pol_dim == settings{i}.dim
               error('Polarization should be parallel, but normal vector directions match.');
            end
            results{i}.gain = (n1*cos(aoi)-n2*sqrt(1-(n1/n2*sin(aoi)).^2)) ./ (n1*cos(aoi)+n2*sqrt(1-(n1/n2*sin(aoi)).^2));
         otherwise
            error('Unsupported partitioning.names "%s".', lower(partitioning.names{j}));
      end
      
      % check mode / modify gain factor (off or transmission); DO NOT USE MODE SETTING 
      % (depends on internalsettings that could be wrong)
      switch lower(results{i}.mode)
         case 'off' % no reflection / transmission 
            results{i}.gain = 0;
         case 'reflection' % do nothing (Fresnel)
         case 'transmission'
            results{i}.gain = sqrt(1 - results{i}.gain.^2).^2; % transmission gain: sqrt(1-refl_gain.^2); air-medium-air: squared
         otherwise
            error('Unsupported results{%i}.mode = %s.', i, results{i}.mode);
      end
      
   end
end


% *******************************************************************************************************
% complete/check partitioning

% remove first entry from partitioning.indices (not needed in test function)
partitioning.indices(1) = [];

% one last check
if ( sum(partitioning.runs) ~= length(results) ) || ( sum(partitioning.runs) ~= length(settings) )
   criterr('Unexpected length of results or settings array. Check partitioning.');
end


% *******************************************************************************************************
% save results

% is there already a .mat-file with this name?
search = dir(strcat(fullfile(directory, filename),'*'));
if ~isempty( search )
   % ask before overwriting
   reply = input('This will overwrite an existing file. Are you REALLY sure y/n [n]?  ', 's');
   if ~strcmpi(reply, 'y')
      return
   end
end

% prepare information header
matfilename    = filename; %#ok<NASGU>
characteristic = 'settings and results for selftest: channel_surface.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
