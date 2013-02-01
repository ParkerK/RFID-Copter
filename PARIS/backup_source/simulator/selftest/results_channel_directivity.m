% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_channel_directivity function
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
filename  = 'results_channel_directivity';
 
% random test setup bounds (uniformly distributed)
testsettings.txant  = {'', 'channelchar_directivity_l4-dipole', 'channelchar_directivity_l2-dipole'}; % transmitter antenna types (empty: uniform)
testsettings.rxant  = {'', 'channelchar_directivity_l4-dipole', 'channelchar_directivity_l2-dipole'}; % receiver antenna types (empty: uniform)
testsettings.ndir   = [1, 10]; % number of gains to be calculated in parallel (size of .dir_tx and .dir_rx)

% partitioning
partitioning.names = {'default: dir_tx', 'option: dir_tx and dir_rx'};
partitioning.runs  = [1000, 1000];


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
      % antenna types
      settings{i}.txant = testsettings.txant{randi(length(testsettings.txant))};
      settings{i}.rxant = testsettings.rxant{randi(length(testsettings.rxant))};
      % antenna rotation [azimuth, elevation]
      %  ... also test the function's ability to deal with nonstd sphere-coordinates (ambiguities)
      settings{i}.txrot = [720*rand-360, 720*rand-360];
      settings{i}.rxrot = [720*rand-360, 720*rand-360];
      % direction from point of transmitter
      ndir = randi(testsettings.ndir);
      %  ... also test the function's ability to deal with nonstd sphere-coordinates (ambiguities)
      %     standard: only dir_tx
      settings{i}.dir_tx = [720*rand(ndir,1)-360, 720*rand(ndir,1)-360];
      %     optional: dir_tx and dir_rx (TX->RX is non-line-of-sight)
      if strcmpi(partitioning.names{j}, 'option: dir_tx and dir_rx')
         settings{i}.dir_rx = [720*rand(ndir,1)-360, 720*rand(ndir,1)-360];
      end
   end
end


% *******************************************************************************************************
% calculate expected results

for i = 1 : length(settings)
   
   % combine path directions and antenna rotations
   %     transmitter
   tx(:,1) = settings{i}.dir_tx(:,1) - settings{i}.txrot(1);
   tx(:,2) = settings{i}.dir_tx(:,2) - settings{i}.txrot(2);
   %     receiver
   if isfield(settings{i}, 'dir_rx')
      rx(:,1) = settings{i}.dir_rx(:,1) - settings{i}.rxrot(1);
      rx(:,2) = settings{i}.dir_rx(:,2) - settings{i}.rxrot(2) - 180;
   else
      rx(:,1) = settings{i}.dir_tx(:,1) - settings{i}.rxrot(1);
      rx(:,2) = settings{i}.dir_tx(:,2) - settings{i}.rxrot(2) - 180;
   end
   
   % get angles to standard sphere coordinates (azimuth 0..2*pi, elevation 0..pi)
   %     transmitter  
   [x,y,z] = sph2cart(tx(:,1)*pi/180, pi/2 - tx(:,2)*pi/180, ones(size(tx,1),1));
   tx(:,1) = mod(atan2(y, x), 2*pi);
   tx(:,2) = pi/2 - atan(z ./ sqrt(x.^2+y.^2));
   %     transmitter
   [x,y,z] = sph2cart(rx(:,1)*pi/180, pi/2 - rx(:,2)*pi/180, ones(size(rx,1),1));
   rx(:,1) = mod(atan2(y, x), 2*pi);
   rx(:,2) = pi/2 - atan(z ./ sqrt(x.^2+y.^2));
   %     check angles
   if any(tx(:,1) < 0) || any(tx(:,1) > 2*pi)
      error('TX azimuth out of range.')
   end
   if any(tx(:,2) < 0) || any(tx(:,2) > pi)
      error('TX elevation out of range.')
   end
   if any(rx(:,1) < 0) || any(rx(:,1) > 2*pi)
      error('RX azimuth out of range.')
   end
   if any(rx(:,2) < 0) || any(rx(:,2) > pi)
      error('RX elevation out of range.')
   end
   %     store angles for later comparison (helps in debugging)
   results{i}.angle_tx = tx * 180/pi;
   results{i}.angle_rx = rx * 180/pi;
   
   % calculate antenna gains
   %     transmitter
   switch lower(settings{i}.txant)
      case '' % uniform
         gain.tx      = ones(size(tx,1), 2);
      case 'channelchar_directivity_l4-dipole' % lambda/4 dipole in z
         gain.tx(:,1) = ones(size(tx(:,1)));
         gain.tx(:,2) = ( cos(pi/2 * cos(tx(:,2))) - cos(pi/2) ) ./ sqrt( 1 - cos(tx(:,2)).^2 );
      case 'channelchar_directivity_l2-dipole' % lambda/2 dipole in z
         gain.tx(:,1) = ones(size(tx(:,1)));
         gain.tx(:,2) = ( cos(pi * cos(tx(:,2))) - cos(pi) ) ./ sqrt( 1 - cos(tx(:,2)).^2 );
      otherwise
         error('Unsupported TX antenna type: %s.', settings{i}.txant);
   end
   %     receiver
   switch lower(settings{i}.rxant)
      case '' % uniform
         gain.rx      = ones(size(rx,1), 2);
      case 'channelchar_directivity_l4-dipole' % lambda/4 dipole in z
         gain.rx(:,1) = ones(size(rx(:,1)));
         gain.rx(:,2) = ( cos(pi/2 * cos(rx(:,2))) - cos(pi/2) ) ./ sqrt( 1 - cos(rx(:,2)).^2 );
      case 'channelchar_directivity_l2-dipole' % lambda/2 dipole in z
         gain.rx(:,1) = ones(size(rx(:,1)));
         gain.rx(:,2) = ( cos(pi * cos(rx(:,2))) - cos(pi) ) ./ sqrt( 1 - cos(rx(:,2)).^2 );
      otherwise
         error('Unsupported RX antenna type: %s.', settings{i}.rxant);
   end
   %     combine
   results{i}.gain = prod(gain.tx,2) .* prod(gain.rx,2);
   
   % cleanup before next setting
   clear('tx','rx','gain');
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
characteristic = 'settings and results for selftest: channel_directivity.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
