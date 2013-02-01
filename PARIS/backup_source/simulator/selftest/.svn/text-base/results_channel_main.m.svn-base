% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_channel_main function
% !!! EXISTING RESULTS WILL BE OVERWRITTEN - USE WITH CARE !!!
%
% Note that the settings for all subfunctions are chosen randomly and independent of each other, so some
% settings might be (geometrically) contradictory!
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

% initialize global stuff (assumes that globalinit.m is part of path)
globalinit('silent');

% "switch to here"
cd(fileparts(mfilename('fullpath')));


% *******************************************************************************************************
% test setup

% output directory and filename
directory = 'results';
filename  = 'results_channel_main';

% nonrandom setup
%     general
stdsettings.noise_mintxlen = 2e4; % samples minimum transmit signal length if noise is added (for estimators)
stdsettings.tol_gain       = eps; % tolerance for gain/delay combination (no identical delay entries) check
%     channel_large
stdsettings.large   = 'log-dist'; % log-distance large-scale model
%     channel_small
stdsettings.maxiter = 1000; % maximum number of iterations for PDP generation
stdsettings.maxrays =  100; % maximum number of paths
%     channel_noise
stdsettings.noise = 'wgn'; % white, Gaussian noise
 
% random test setup bounds (uniformly distributed)
%     general settings
testsettings.fs    = [  1e9, 10e9]; % Hz sampling frequency
testsettings.txa   = [  0.1,   10]; % sent signal amplitudes
testsettings.ntx   = [    1,    3]; % number of transmitters
testsettings.nrx   = [    1,    3]; % number of receivers
%     channel_main
testsettings.p_dls =         0.1; % probability for direct large-scale gain/delay
testsettings.dls_g = [1e-4,   1]; % direct large-scale gain fayctor
testsettings.dls_d = [   1, 1e3]; % direct large-scale delay (in samples)
%     channel_large
testsettings.dist = [  1,   10]; % m distance
testsettings.pl   = [1.5,    4]; % path-loss factor
testsettings.c    = [3e7,  8e7]; % m/s speed of light ("tired light" to get higher delay = longer vectors)
testsettings.f0   = [1e6, 30e6]; % Hz "center" frequencies for attenuation
%     channel_directivity
testsettings.txant  = {'', 'channelchar_directivity_l4-dipole', 'channelchar_directivity_l2-dipole'}; % transmitter antenna types
testsettings.rxant  = {'', 'channelchar_directivity_l4-dipole', 'channelchar_directivity_l2-dipole'}; % receiver antenna types
%     channel_reflection
testsettings.surf = {'floor', 'wall', 'ceiling'}; % a few surfaces
testsettings.n2   =                      [1, 10]; % refractive index of surface
%     channel_small (only deterministic setup => seed has no effect)
testsettings.k        = [   10,   30]; % dB Ricean K-factor
testsettings.trms     = [ 2e-9, 2e-8]; % s RMS delay spread
testsettings.bw       = [  1e6,  1e8]; % Hz bandwidth for channel (with oversampling)
testsettings.eps_k    = [ 5e-4, 2e-3]; % relative tolerable error in K-factor (linear)
testsettings.eps_trms = [ 5e-3, 2e-2]; % relative tolerable error in RMS delay spread
testsettings.fres     = [ 1e-7, 5e-7]; % Hz initial frequency resolution
%     channel_noise
testsettings.n0     = [-80, -40; 1, 1e3]; % [dBm; Hz] single-sided noise density [dBm, per ? Hz] MEASURED AT RECEIVER
testsettings.frxs   = [1e6, 10e6]; % Hz receiver (tag,reader) sampling frequency (for noise density -> power)

% partitioning
partitioning.names = {...
   'largescale/direct', 'largescale/direct + noise', 'largescale/direct + reflections',...
   'largescale/direct + smallscale', 'largescale/direct + reflections + smallscale',...
   'largescale/direct + reflections + smallscale + noise'};
partitioning.runs  = [500, 200, 300, 100, 100, 20];


% *******************************************************************************************************
% complete/check partitioning

if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check partitioning.');
end

% indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% check channel_large, channel_directivity, channel_reflection, and channel_small (used below)

disp('Starting test routines to make sure channel_large, channel_directivity, channel_reflection, and channel_small are functional.');

errors = test_channel_large(0);
errors = test_channel_directivity(errors);
errors = test_channel_reflection(errors);
errors = test_channel_small(errors);
if errors ~= 0
   disp(sprintf('\n*** TEST FAILED, ABORTING ***'));
   return
end


% *******************************************************************************************************
% create settings
disp(sprintf('\nCreating settings...'));

% noiseless
for i = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)
      
      % general
      %     sampling frequency
      settings{j}.fs = testsettings.fs(1) + rand *(testsettings.fs(2) - testsettings.fs(1));
      %     number of transmitters, receivers
      settings{j}.ntx = randi(testsettings.ntx);
      settings{j}.nrx = randi(testsettings.nrx);
      %     switches
      settings{j}.switches.direct_gain = round(rand-0.5+testsettings.p_dls);
      switch lower(partitioning.names{i})
         case 'largescale/direct'
            settings{j}.switches.refl_on   = false; % reflective surfaces
            settings{j}.switches.small_on  = false; % smallscale channel
            settings{j}.switches.noise_on  = false; % noise
         case 'largescale/direct + noise'
            settings{j}.switches.refl_on   = false;
            settings{j}.switches.small_on  = false;
            settings{j}.switches.noise_on  =  true;
         case 'largescale/direct + reflections'
            settings{j}.switches.refl_on   =  true;
            settings{j}.switches.small_on  = false;
            settings{j}.switches.noise_on  = false;
         case 'largescale/direct + smallscale'
            settings{j}.switches.refl_on   = false;
            settings{j}.switches.small_on  =  true;
            settings{j}.switches.noise_on  = false;
         case 'largescale/direct + reflections + smallscale'
            settings{j}.switches.refl_on   =  true;
            settings{j}.switches.small_on  =  true;
            settings{j}.switches.noise_on  = false;
         case 'largescale/direct + reflections + smallscale + noise'
            settings{j}.switches.refl_on   =  true;
            settings{j}.switches.small_on  =  true;
            settings{j}.switches.noise_on  =  true;
         otherwise
            error('Unsupported partitioning.names "%s".', lower(partitioning.names{i}));
      end
      %     reflective surfaces
      if settings{j}.switches.refl_on
         num_refl = randi(length(testsettings.surf));
         settings{j}.ind_refl = randperm(length(testsettings.surf));
         settings{j}.ind_refl = settings{j}.ind_refl(1:num_refl);
      end
      clear('num_refl');

      % for transmitters
      for tx = 1 : settings{j}. ntx
         % sent signal amplitudes
         settings{j}.channel{tx,1}.txa = round(testsettings.txa(1) + rand*(testsettings.txa(2) - testsettings.txa(1)));
         
         % for receivers
         for rx = 1 : settings{j}.nrx
            % copy sent signal amplitudes to all other receivers (just to make sure)
            settings{j}.channel{tx,rx}.txa = settings{j}.channel{tx,1}.txa;

            % large-scale and directivity
            %     direct gain/delay
            if settings{j}.switches.direct_gain
               settings{j}.channel{tx,rx}.direct.gain    = ...
                        testsettings.dls_g(1) + rand*(testsettings.dls_g(2) - testsettings.dls_g(1) );
               settings{j}.channel{tx,rx}.direct.delay_s = ...
                  round(testsettings.dls_d(1) + rand*(testsettings.dls_d(2) - testsettings.dls_d(1) ));

            %     normal largescale model
            else
               % distance matrix (geometry neglected)
               settings{j}.channel{tx,rx}.largescale.dist = testsettings.dist(1) + rand*(testsettings.dist(2) - testsettings.dist(1) );
               % other settings
               settings{j}.channel{tx,rx}.largescale.type = stdsettings.large;
               settings{j}.channel{tx,rx}.largescale.pl   = testsettings.pl(1) + rand*(testsettings.pl(2) - testsettings.pl(1));
               settings{j}.channel{tx,rx}.largescale.c    = testsettings.c(1)  + rand*(testsettings.c(2)  - testsettings.c(1) );
               settings{j}.channel{tx,rx}.largescale.f0   = testsettings.f0(1) + rand*(testsettings.f0(2) - testsettings.f0(1));
               settings{j}.channel{tx,rx}.largescale.fs   =  settings{j}.fs;
               % round distances to sampling resolution
               settings{j}.channel{tx,rx}.largescale.dist = round(settings{j}.channel{tx,rx}.largescale.dist/...
                  settings{j}.channel{tx,rx}.largescale.c*settings{j}.channel{tx,rx}.largescale.fs)*...
                  settings{j}.channel{tx,rx}.largescale.c/settings{j}.channel{tx,rx}.largescale.fs;
               
               % directivity (antenna gain patterns)
               settings{j}.channel{tx,rx}.directivity.txant  = testsettings.txant{randi(length(testsettings.txant))};
               settings{j}.channel{tx,rx}.directivity.rxant  = testsettings.rxant{randi(length(testsettings.rxant))};
               settings{j}.channel{tx,rx}.directivity.txrot  = [360*rand, 180*rand];
               settings{j}.channel{tx,rx}.directivity.rxrot  = [360*rand, 180*rand];
               settings{j}.channel{tx,rx}.directivity.dir_tx = [360*rand, 180*rand];
            end
            
            % reflective surfaces
            if settings{j}.switches.refl_on
               for k = 1 : length(settings{j}.ind_refl)
                  refl_name = testsettings.surf{settings{j}.ind_refl(k)};
                  % reflection coefficient
                  settings{j}.channel{tx,rx}.reflections.(refl_name).aoi     = 90 * rand;
                  settings{j}.channel{tx,rx}.reflections.(refl_name).n2      = testsettings.n2(1) + diff(testsettings.n2) * rand;
                  settings{j}.channel{tx,rx}.reflections.(refl_name).dim     = randi([1,3]); % "3D"
                  settings{j}.channel{tx,rx}.reflections.(refl_name).pol_dim = randi([1,3]); % "3D"
                  % another set of largescale/directivity settings (copied from above)
                  %     largescale
                  settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.dist = ...
                     testsettings.dist(1) + rand*(testsettings.dist(2) - testsettings.dist(1) );
                  settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.type = stdsettings.large;
                  settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.pl   = ...
                     testsettings.pl(1) + rand*(testsettings.pl(2) - testsettings.pl(1));
                  settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.c    = ...
                     testsettings.c(1)  + rand*(testsettings.c(2)  - testsettings.c(1) );
                  settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.f0   = ...
                     testsettings.f0(1) + rand*(testsettings.f0(2) - testsettings.f0(1));
                  settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.fs   =  settings{j}.fs;
                  settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.dist = ...
                     round(settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.dist / ...
                     settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.c * ...
                     settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.fs) * ...
                     settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.c / ...
                     settings{j}.channel{tx,rx}.reflections.(refl_name).largescale.fs;
                  %     directivity
                  settings{j}.channel{tx,rx}.reflections.(refl_name).directivity.txant  = ...
                     testsettings.txant{randi(length(testsettings.txant))};
                  settings{j}.channel{tx,rx}.reflections.(refl_name).directivity.rxant  = ...
                     testsettings.rxant{randi(length(testsettings.rxant))};
                  settings{j}.channel{tx,rx}.reflections.(refl_name).directivity.txrot  = [360*rand, 180*rand];
                  settings{j}.channel{tx,rx}.reflections.(refl_name).directivity.rxrot  = [360*rand, 180*rand];
                  settings{j}.channel{tx,rx}.reflections.(refl_name).directivity.dir_tx = [360*rand, 180*rand];
               end
               % check
               if length(fieldnames(settings{j}.channel{tx,rx}.reflections)) ~= length(settings{j}.ind_refl)
                  error('We seem to have lost a few reflective surfaces...');
               end
            end       

            % small-scale
            %     on/off, number of paths
            settings{j}.channel{tx,rx}.smallscale.on        = settings{j}.switches.small_on;
            settings{j}.channel{tx,rx}.smallscale.det       = true; % only deterministic
            settings{j}.channel{tx,rx}.smallscale.ensembles =    1; % only deterministic
            settings{j}.channel{tx,rx}.smallscale.seed      =  NaN; % only deterministic => has no effect anyway
            settings{j}.channel{tx,rx}.smallscale.maxrays   = stdsettings.maxrays; % with given maximum number of paths
            settings{j}.channel{tx,rx}.smallscale.maxiter   = stdsettings.maxiter; % and iterations
            %     other settings
            settings{j}.channel{tx,rx}.smallscale.k        = testsettings.k(1)        + rand*(testsettings.k(2)        - testsettings.k(1)       );
            settings{j}.channel{tx,rx}.smallscale.trms     = testsettings.trms(1)     + rand*(testsettings.trms(2)     - testsettings.trms(1)    );
            settings{j}.channel{tx,rx}.smallscale.bw       = testsettings.bw(1)       + rand*(testsettings.bw(2)       - testsettings.bw(1)      );
            settings{j}.channel{tx,rx}.smallscale.eps_k    = testsettings.eps_k(1)    + rand*(testsettings.eps_k(2)    - testsettings.eps_k(1)   );
            settings{j}.channel{tx,rx}.smallscale.eps_trms = testsettings.eps_trms(1) + rand*(testsettings.eps_trms(2) - testsettings.eps_trms(1));
            settings{j}.channel{tx,rx}.smallscale.fres     = testsettings.fres(1)     + rand*(testsettings.fres(2)     - testsettings.fres(1)    );
            settings{j}.channel{tx,rx}.smallscale.fs       =  settings{j}.fs;
         end
      end
      
      % noise: for receivers
      for rx = 1 : settings{j}.nrx
         % set first transmitter randomly
         %     on
         if settings{j}.switches.noise_on
            settings{j}.channel{1,rx}.noise.type = stdsettings.noise;
            settings{j}.channel{1,rx}.noise.fs   = settings{j}.fs;
            settings{j}.channel{1,rx}.noise.n0   = [testsettings.n0(1,1) + rand*(testsettings.n0(1,2) - testsettings.n0(1,1)),...
                                                    testsettings.n0(2,1) + rand*(testsettings.n0(2,2) - testsettings.n0(2,1))];
            settings{j}.channel{1,rx}.noise.frxs =  testsettings.frxs(1) + rand*(testsettings.frxs(2) - testsettings.frxs(1));
            %     off (settings have to be complete nontheless)
         else
            settings{j}.channel{1,rx}.noise.type = 'off';
            settings{j}.channel{1,rx}.noise.n0   = [-Inf,1]; % -Inf dB per Hz => no noise
            settings{j}.channel{1,rx}.noise.fs   = 1;
            settings{j}.channel{1,rx}.noise.frxs = 1;
         end

         % copy all others (just to make sure)
         for tx = 1 : settings{j}.ntx
            settings{j}.channel{tx,rx}.noise.type = settings{j}.channel{1,rx}.noise.type;
            settings{j}.channel{tx,rx}.noise.n0   = settings{j}.channel{1,rx}.noise.n0;
            settings{j}.channel{tx,rx}.noise.fs   = settings{j}.channel{1,rx}.noise.fs;
            settings{j}.channel{tx,rx}.noise.frxs = settings{j}.channel{1,rx}.noise.frxs;
         end
      end
   end
end

% clean up
clear('i', 'j', 'k', 'num_refl', 'refl_name', 'rx', 'tx');


% *******************************************************************************************************
% calculate expected results and set input signal lengths (minimum possible)
disp(sprintf('Creating results...'));

% calculate results
for i = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)
      
      % generate filter channels
      for tx = 1 : settings{j}.ntx
         for rx = 1 : settings{j}.nrx
            
            % main channel model
            %     largescale (wave propagation)
            if settings{j}.switches.direct_gain % direct manipulation of gain/delay
               delay_ls = settings{j}.channel{tx,rx}.direct.delay_s;
               gain_ls  = settings{j}.channel{tx,rx}.direct.gain;
            else % normal largescale model
               [delay_ls, gain_ls] = channel_large(settings{j}.channel{tx,rx}.largescale);
               gain_ls = gain_ls * channel_directivity(settings{j}.channel{tx,rx}.directivity);
            end
            %     smallscale (statistical model)
            if settings{j}.switches.small_on
               [delay_ss, gain_ss] = channel_small(settings{j}.channel{tx,rx}.smallscale);
            end
            
            % reflective surfaces
            if settings{j}.switches.refl_on
               delay_rs = nan(length(settings{j}.ind_refl), 1);
               gain_rs = nan(length(settings{j}.ind_refl), 1);
               for k = 1 : length(settings{j}.ind_refl)
                  refl_name = testsettings.surf{settings{j}.ind_refl(k)};
                  % wave propagation
                  [delay_rs(k), gain_rs(k)] = ...
                     channel_large(settings{j}.channel{tx,rx}.reflections.(refl_name).largescale);
                  % antenna gain (directivity)
                  gain_rs(k) = gain_rs(k) * ...
                     channel_directivity(settings{j}.channel{tx,rx}.reflections.(refl_name).directivity);
                  % reflection coefficient of reflective surface
                  gain_rs(k) = gain_rs(k) * ...
                     channel_reflection(settings{j}.channel{tx,rx}.reflections.(refl_name));
               end
            end
            
            % combine filter channels
            %     main channel (largescale + smallscale)
            results{j}.delay_s{tx,rx} = delay_ls;
            results{j}.gain{tx,rx} = gain_ls;
            if settings{j}.switches.small_on
               results{j}.delay_s{tx,rx} = results{j}.delay_s{tx,rx} + delay_ss;
               results{j}.gain{tx,rx} = results{j}.gain{tx,rx} * gain_ss;
            end
            %     deterministic reflections
            if settings{j}.switches.refl_on
               results{j}.delay_s{tx,rx} = [results{j}.delay_s{tx,rx}; delay_rs];
               results{j}.gain{tx,rx} =    [results{j}.gain{tx,rx}; gain_rs];
               % sort and add up identical entries
               [results{j}.delay_s{tx,rx}, sort_ind] = sort(results{j}.delay_s{tx,rx});
               results{j}.gain{tx,rx} = results{j}.gain{tx,rx}(sort_ind);
               %     construct a delay vector without identical entries
               clean_delay = results{j}.delay_s{tx,rx};
               clean_delay(2:end) = diff(clean_delay);
               clean_delay([false; clean_delay(2:end)==0]) = []; % make sure first entry is not deleted
               clean_delay = cumsum(clean_delay);
               %     from that, construct a clean gain vector
               clean_gain = nan(size(clean_delay));
               for k = 1 : length(clean_delay)
                  clean_gain(k) = sum( results{j}.gain{tx,rx}(results{j}.delay_s{tx,rx} == clean_delay(k)) );
               end
               %     double-check result
               if sum(results{j}.gain{tx,rx}) - sum(clean_gain) > stdsettings.tol_gain
                  error('Sum of gains does not match (after simplification of filter impulse response)');
               end
               %     and replace originals
               results{j}.gain{tx,rx} = clean_gain;
               results{j}.delay_s{tx,rx} = clean_delay;
            end
            
         end
         
         % test vector lengths (maximum delay; delay_s is sorted)
         settings{j}.tx_len(tx) = max(cellfun(@(x) x(end), results{j}.delay_s(tx,:))) + 1; % delay=0 for small only
         %     if noise: minimum length for estimators
         if settings{j}.switches.noise_on
            settings{j}.tx_len(tx) = max(settings{j}.tx_len(tx), stdsettings.noise_mintxlen);
         end
      end
      
      % all RX signals should have length identical to maximum input signal length; get maximum length
      results{j}.rx_len = max(settings{j}.tx_len);
        
      % create overlay (delay_s == 0 => index == 1)
      %     get maximum input signal length and create matrix (all output lengths are identical)
      overlay = zeros(results{j}.rx_len, settings{j}.nrx);
      %     vectorize delays and gains and add up gain factors (does not have to be performant...)
      for rx = 1 : settings{j}.nrx % for each receiver
         for tx =  1 : settings{j}.ntx % for each transmitter
            overlay(results{j}.delay_s{tx,rx}+1, rx) = overlay(results{j}.delay_s{tx,rx}+1, rx) +...
               squeeze(results{j}.gain{tx,rx}) * settings{j}.channel{tx,rx}.txa;
         end
      end
      %     this will likely be a quite sparse matrix...
      results{j}.ol = sparse(overlay);
      
      % calculate noise variances for receivers (settings should be identical for all transmitters)
      for rx = 1 : settings{j}.nrx
         results{j}.nvar(rx) = 10^(settings{j}.channel{1,rx}.noise.n0(1)/10-3)/settings{j}.channel{1,rx}.noise.n0(2) *...
            settings{j}.channel{1,rx}.noise.fs./settings{j}.channel{1,rx}.noise.frxs; % n0: [dBm; Hz]
      end
      
   end
end


% *******************************************************************************************************
% complete/check partitioning

% remove first entry from partitioning.indices (not needed in test function)
partitioning.indices(1) = [];

% one last check
if ( sum(partitioning.runs) ~= length(results) ) || ( sum(partitioning.runs) ~= length(settings) )
   criterr('Unexpected length of settings array. Check partitioning.');
end


% *******************************************************************************************************
% save results

% empty line
disp(sprintf('\n'));

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
characteristic = 'settings and results for selftest: channel_main.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));

