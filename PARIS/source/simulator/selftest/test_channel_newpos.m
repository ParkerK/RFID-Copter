% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% selftest: channel - adjust all channel parameters to new transmitter/receiver positions
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
% version = test_channel_newpos()
%    Just returns the version number (string).
% sumoferrors = test_channel_newpos(sumoferrors)
%    Tests the function CHANNEL_NEWPOS, displays the results in the command window and returns the 
%    overall amount of errors found.
%
%
% ***** Interface definition *****
% function sumoferrors = test_channel_newpos(sumoferrors);
%    sumoferrors    sum of errors found not including this fcn call
%
%    sumoferrors    sum of errors found including this fcn call (+1 for each erroneous tested
%                   functionality)
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
% = re-implement after re-implementation of channel_newpos
%
% *******************************************************************************************************


function sumoferrors = test_channel_newpos(sumoferrors)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   sumoferrors = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% tolerances
%     distance
internalsettings.tol_dist  = 1e-4; % relative error
%     angles (may be zero)
internalsettings.tol_angle = 1e-4; % degree
%     K-factor and RMS delay spread
internalsettings.tol_sspar = 1e-4; % relative error


% *******************************************************************************************************
% initialize

% output
disp('   = channel_newpos ***');
error('     TEST IS OUTDATED');
disp(sprintf('     checking distances +/-%g%%, angles +/-%gdeg',...
   internalsettings.tol_dist*100, internalsettings.tol_angle));

% load settings and expected results
data = loadmat('results_channel_newpos.mat', '     ');
settings     = data.settings;
partitioning = data.partitioning;
results      = data.results;


% *******************************************************************************************************
% run tests and evaluate

for i = 1 : length(results)
   
   % partitioning and output
   %     find index for partition struct (partitioning.indices is "end-of-block")
   for j = 1 : length(partitioning.indices)
      if i <= partitioning.indices(j)
         index = j;
         break;
      end
   end
   
   % display once / reset errors for each new test block
   if i == 1 || ( index > 1 && i == partitioning.indices(index-1) + 1 )
      disp(sprintf('      - %i setup(s) "%s"', partitioning.runs(index), partitioning.names{index} ));
      if strcmpi(partitioning.names{index}, 'los/direct')
         disp(sprintf('        additional check: K-factor and RMS delay spread +/-%g%%', internalsettings.tol_sspar*100));
      end
      % reset errors
      errors    = 0;
      errortext = '';
   end
   
   % create call parameters for channel_newpos
   %     channel parameters to modify
   settings2mod = cell(settings{i}.ntx, settings{i}.nrx);
   for tx = 1 : settings{i}.ntx
      for rx = 1 : settings{i}.nrx
         % largescale/directivity or feedback channel?
         if settings{i}.no_ls(tx,rx)
            settings2mod{tx,rx} = struct();
         else
            settings2mod{tx,rx} = struct('largescale',struct(), 'directivity',struct());
         end
         % virtual connections?
         if settings{i}.is_vtx(tx,rx)
            settings2mod{tx,rx}.virtual.vtx = true;
            settings2mod{tx,rx}.virtual.otx = settings{i}.otx(tx,rx);
         end
         if settings{i}.is_vrx(tx,rx)
            settings2mod{tx,rx}.virtual.vrx = true;
            settings2mod{tx,rx}.virtual.orx = settings{i}.orx(tx,rx);
         end
      end
   end
   %     settings for channel_newpos
   ownsettings = struct('v_dist',settings{i}.v_dist, 'v_k',settings{i}.v_k, 'v_trms',settings{i}.v_trms);
   if isfield(settings{i}, 'reflections')
      ownsettings.reflections = settings{i}.reflections;
   end
   
   % call channel_newpos  
   settings2mod = channel_newpos(settings2mod, settings{i}.pos_t, settings{i}.pos_r, ownsettings);
   
   % check presence of necessary fields
   testresults.no_largescale  = cellfun(@(x) ~isfield(x, 'largescale'), settings2mod);
   testresults.no_directivity = cellfun(@(x) ~isfield(x, 'directivity'), settings2mod);
   [errors, errortext] = check_result(i, testresults.no_largescale,  settings{i}.no_ls, ...
      errors, errortext, 'largescale missing', 'equal'); 
   [errors, errortext] = check_result(i, testresults.no_directivity, settings{i}.no_ls, ...
      errors, errortext, 'directivity missing', 'equal');
   
   % extract relevant data (random checks in order to keep results_channel_newpos as simple as possible)
   testresults.dist   = nan(settings{i}.nrx, 1);
   testresults.dir_tx = nan(settings{i}.nrx, 2);
   testresults.k      = nan(settings{i}.nrx, 1);
   testresults.trms   = nan(settings{i}.nrx, 1);
   for rx = 1 : settings{i}.nrx
      % limit to one transmitter and one reflective surface
      if isfield(settings{i}, 'reflections')
         temp = settings2mod{results{i}.ind_t0, rx}.reflections.(settings{i}.refl{settings{i}.ind_refl0});
      else
         temp = settings2mod{results{i}.ind_t0, rx};
      end
      if isfield(temp, 'largescale')
         testresults.dist(rx) = temp.largescale.dist;
      end
      if isfield(temp, 'directivity')
         testresults.dir_tx(rx,:) = mod(temp.directivity.dir_tx, 360);
      end
      testresults.k(rx)    = settings2mod{results{i}.ind_t0, rx}.smallscale.k;
      testresults.trms(rx) = settings2mod{results{i}.ind_t0, rx}.smallscale.trms;
   end
   
   % checks
   %     largescale/directivity (also for reflection)
   [errors, errortext] = check_result(i, testresults.dist,   results{i}.dir_t0r(:,end),...
      errors, errortext, 'dist', 'relerr', internalsettings.tol_dist); 
   [errors, errortext] = check_result(i, testresults.dir_tx(:,1), results{i}.dir_t0r(:,1),...
      errors, errortext, 'dir_tx(az)', 'range', [-1,1]*internalsettings.tol_angle);
   [errors, errortext] = check_result(i, testresults.dir_tx(:,2), results{i}.dir_t0r(:,2),...
      errors, errortext, 'dir_tx(el)', 'range', [-1,1]*internalsettings.tol_angle);
   %     smallscale (depend on LOS distance => invalid if we're dealing with a reflection)
   if ~isfield(settings{i}, 'reflections')
      [errors, errortext] = check_result(i, testresults.k,results{i}.k,...
         errors, errortext, 'K',    'relerr', internalsettings.tol_sspar);
      [errors, errortext] = check_result(i, testresults.trms, results{i}.trms,...
         errors, errortext, 'trms', 'relerr', internalsettings.tol_sspar);
   end
         
   % output (if current loop is last test in block)
   if i == partitioning.indices(index)
      if errors == 0
         disp('         ... passed');
      else
         disp(sprintf('         ... ERRORS: %s', errortext));
         sumoferrors = sumoferrors + 1;  
      end
   end   
end

