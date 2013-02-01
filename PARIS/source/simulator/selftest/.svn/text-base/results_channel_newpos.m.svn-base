% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_channel_newpos function
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
% = massively outdated; re-implement after re-implementation of channel_newpos; 
%   add missing parameters (AOI, surfaces, ...)
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
filename  = 'results_channel_newpos';

% nonrandom test setup
%     smallscale parameters
testsettings.npt_ss = 2; % number of points for parameter vectors (can be small if poly* below are linear)
%     ... polynomials: linear for best results w. linear interp. (make sure K, trms >> 0)
testsettings.maxdist =       1000; % m assumed maximum distance between objects (for smallscale parameters)
testsettings.poly_k  = [0.31 222]; % polynomial for K-factor vs. distance (see polyval)
testsettings.poly_t  = [0.15  33]; % polynomial for RMS delay spread vector vs. distance (see polyval) 

% random test setup bounds [min, max] (uniformly distributed)
%     number of...
testsettings.ndim = [  1,  3]; % number of dimensions (keep at 1..3)
testsettings.ntx  = [  1, 10]; % number of transmitters (choose max. < nrx, only nrx connections checked!)
testsettings.nrx  = [  1, 30]; % number of receivers
%     placement, distances
testsettings.xyz  = [-10, 10]; % m bounds for placement of first transmitter/receiver (all dimensions)
testsettings.dist = [0.1, 10]; % m maximum distance for each step of the positions trajectory
%     reflective surfaces
testsettings.surf  = {'floor', 'wall', 'ceiling'}; % a few surfaces
testsettings.shift = [-10, 10]; % m maximum shift of reflective surface along .dim (surface normal)
%     probability that there is no largescale/directivity model (feedback channel)
testsettings.p_nols = 0.05;
%     probability that this is a VTX or VRX connection (VTX, VRX equally distributed)
testsettings.p_virt = 0.75;        

% partitioning
partitioning.names = {'Direct', 'Reflections', 'Reflections+Virtual'};
partitioning.runs  = [200, 100, 100];

% internal checks
%     for testsettings.p_nlos
internal.bds_nols = testsettings.p_nols/10; % maximum deviation to testsettings.p_nols
internal.n_ls     =                      0; % number of normal channels
internal.n_nols   =                      0; % number of direct feedback channels


% *******************************************************************************************************
% complete/check partitioning

if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check partitioning.');
end

% indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% create settings / (some) expected results

for j = 1 : length(partitioning.names)
   for i = 1 + partitioning.indices(j) : partitioning.indices(j+1)
      % number of...
      settings{i}.ndim = randi(testsettings.ndim); % dimensions
      settings{i}.ntx  = randi(testsettings.ntx);  % transmitters
      settings{i}.nrx  = randi(testsettings.nrx);  % receivers
      %     abbreviations for below
      ndim = settings{i}.ndim; 
      ntx  = settings{i}.ntx; 
      nrx  = settings{i}.nrx;

      % positions
      %     random transmitter positions
      settings{i}.pos_t = testsettings.xyz(1) + diff(testsettings.xyz) * rand(ntx, ndim); %#ok<*SAGROW>
      %     random distances from one random transmitter to all receivers (simpler to check)
      results{i}.ind_t0 = randi(ntx);
      switch ndim
         case 1
            results{i}.dir_t0r  = [zeros(nrx,1), pi/2*ones(nrx,1),...
               testsettings.dist(1) + diff(testsettings.dist) * rand(nrx,1)]; % [az, el, dist]
            [dx,dy,dz] = sph2cart(results{i}.dir_t0r(:,1), pi/2-results{i}.dir_t0r(:,2), results{i}.dir_t0r(:,3));
            settings{i}.pos_r = repmat(settings{i}.pos_t(results{i}.ind_t0,:), nrx, 1) + dx;
         case 2
            results{i}.dir_t0r  = [2*pi*rand(nrx,1), pi/2*ones(nrx,1),...
               testsettings.dist(1) + diff(testsettings.dist) * rand(nrx,1)]; % [az, el, dist]
            [dx,dy,dz] = sph2cart(results{i}.dir_t0r(:,1), pi/2-results{i}.dir_t0r(:,2), results{i}.dir_t0r(:,3));
            settings{i}.pos_r = repmat(settings{i}.pos_t(results{i}.ind_t0,:), nrx, 1) + [dx, dy];
         case 3
            results{i}.dir_t0r  = [2*pi*rand(nrx,1), pi*rand(nrx,1),...
               testsettings.dist(1) + diff(testsettings.dist) * rand(nrx,1)]; % [az, el, dist]
            [dx,dy,dz] = sph2cart(results{i}.dir_t0r(:,1), pi/2-results{i}.dir_t0r(:,2), results{i}.dir_t0r(:,3));
            settings{i}.pos_r = repmat(settings{i}.pos_t(results{i}.ind_t0,:), nrx, 1) + [dx, dy, dz];
         otherwise
            error('This number of dimensions (%i) is not supported.', ndim)
      end
      %     for reflective surfaces
      %       ... note that only the results for one out of all surface will be checked
      if strcmpi(partitioning.names{j}, 'Reflections') || strcmpi(partitioning.names{j}, 'Reflections+Virtual')
         % define some reflective surfaces
         num_refl = randi(length(testsettings.surf));
         ind_refl = randperm(length(testsettings.surf));
         ind_refl = ind_refl(1:num_refl);
         % settings for all surfaces
         for k = 1 : length(ind_refl)
            settings{i}.refl{k} = testsettings.surf{ind_refl(k)};
            settings{i}.reflections.(settings{i}.refl{k}).dim   = randi([1,ndim]);
            settings{i}.reflections.(settings{i}.refl{k}).shift = ...
               testsettings.shift(1) + diff(testsettings.shift) * rand;
         end
         % select one and mirror receiver positions => results{i}.dir_t0r valid for this surface
         settings{i}.ind_refl0 = randi([1,length(ind_refl)]);
         settings{i}.pos_r(:, settings{i}.reflections.(settings{i}.refl{settings{i}.ind_refl0}).dim) = ...
            2 * settings{i}.reflections.(settings{i}.refl{settings{i}.ind_refl0}).shift - ...
            settings{i}.pos_r(:, settings{i}.reflections.(settings{i}.refl{settings{i}.ind_refl0}).dim);
      end
      %     channel_newpos is supposed to work for cell arrays and matrices
      if round(rand)
         settings{i}.pos_t = mat2cell(settings{i}.pos_t, ones(ntx,1), ndim);
      end
      if round(rand)
         settings{i}.pos_r = mat2cell(settings{i}.pos_r, ones(nrx,1), ndim);
      end
      
      % smallscale channel parameter vectors
      settings{i}.v_dist = linspace(0, testsettings.maxdist, testsettings.npt_ss);
      settings{i}.v_k    = polyval(testsettings.poly_k, settings{i}.v_dist);
      settings{i}.v_trms = polyval(testsettings.poly_t, settings{i}.v_dist);
      
      % miscellanous
      %     no largescale/directivity model?
      switch lower(partitioning.names{j})
         case 'direct' % no largescale model => direct feedback
            settings{i}.no_ls = rand(ntx, nrx) < testsettings.p_nols;
            settings{i}.is_vtx = false(ntx, nrx);
            settings{i}.is_vrx = false(ntx, nrx);
            %     for internal checks
            internal.n_nols = internal.n_nols + sum(sum(settings{i}.no_ls));
            internal.n_ls   = internal.n_ls + ntx * nrx - sum(sum(settings{i}.no_ls));
         case 'reflections' % NLOS implies that there is a largescale model!
            settings{i}.no_ls  = false(ntx, nrx);
            settings{i}.is_vtx = false(ntx, nrx);
            settings{i}.is_vrx = false(ntx, nrx);
         case 'reflections+virtual'
            settings{i}.no_ls  = false(ntx, nrx);
            is_virt = rand(ntx, nrx) < testsettings.p_virt;
            settings{i}.is_vtx = is_virt & rand(ntx, nrx) < 0.5;
            settings{i}.is_vrx = is_virt & ~settings{i}.is_vtx; % VTX+VRX not supported
            if ~all( (settings{i}.is_vtx(:) | settings{i}.is_vrx(:)) == is_virt(:) )
               error('Problems assigning virtual transmitters / receivers.');
            end
            settings{i}.otx = randi([1,ntx], [ntx,nrx]); % originating transmitter (for VTX)
            settings{i}.otx(~settings{i}.is_vtx) = NaN;
            settings{i}.orx = randi([1,nrx], [ntx,nrx]); % originating receiver (for VRX)
            settings{i}.orx(~settings{i}.is_vrx) = NaN;
         otherwise
            error('Unsupported partitioning.names "%s".', lower(partitioning.names{j}));
      end
      

   end
end

% internal checks
%     testsettings.p_nols
internal.p_nols = internal.n_nols / (internal.n_nols + internal.n_ls);
if abs(internal.p_nols - testsettings.p_nols) > internal.bds_nols
   disp(sprintf('WARNING: testsettings.p_nols was not met (%g vs %g). Small number of tests?', ...
      testsettings.p_nols, internal.p_nols));
end


% *******************************************************************************************************
% calculate missing expected results

for j = 1 : length(partitioning.names)
   for i = 1 + partitioning.indices(j) : partitioning.indices(j+1)
      % direction: angles rad -> deg
      results{i}.dir_t0r(:,1:end-1) = results{i}.dir_t0r(:,1:end-1) * 180/pi;
      %     set largescale/directivity results to NaN if the models are not supposed to exist
      results{i}.dir_t0r(settings{i}.no_ls(results{i}.ind_t0,:), :) = NaN;
      
      % smallscale channel parameters (if supposed to exist)
      switch lower(partitioning.names{j})
         case {'reflections','reflections+virtual'}
            results{i}.k    = polyval(testsettings.poly_k, results{i}.dir_t0r(:,end));
            results{i}.trms = polyval(testsettings.poly_t, results{i}.dir_t0r(:,end));
            %     modifications for virtual transmitters and receivers
            for tx = results{i}.ind_t0 % [TODO: remove this limitation in revised version]
               for rx = 1 : settings{i}.nrx
                  if settings{i}.is_vtx(tx,rx)
                     results{i}.trms(tx,rx) = polyval(testsettings.poly_t, results{i}.dir_t0r(rx,3)); % [all distances calculated from one TX]
                  end
                  if settings{i}.is_vrx(tx,rx)
                     results{i}.trms(tx,rx) = polyval(testsettings.poly_t, results{i}.dir_t0r(settings{i}.orx(tx,rx),3));
                  end
               end
            end
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
characteristic = 'settings and results for selftest: channel_newpos.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
