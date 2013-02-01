% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_mfcw function
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
% ? complete tests such that all parameters are checked (gammai, c_ord, ...)
%
% *******************************************************************************************************


% *******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% initialize global stuff (assumes that globalinit.m is part of path)
globalinit('silent'); % do not print any messages (some functions will produce warnings here)

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% test function intrand()
disp('Running some quick checks on integer random number generator (INTRAND).');
if ~intrand_isok
   error('Integer random number generator does not seem to be functional. Rerun to verify.');
else
   disp('   ...passed.');
end


% *******************************************************************************************************
% test setup

% output directory and filename
directory = 'results';
filename  = 'results_mfcw';

% partitioning
partitioning.names = {'ideal carrier modulation'};
partitioning.runs  = [25];

% standard (nonrandom) setup
%     general setup
stdsettings.fm_per = 48; % number of periods of modulation 
%     channel: largescale
stdsettings.channel{1,1}.largescale.type  = 'log-dist'; % {'log-dist'}
stdsettings.channel{1,1}.largescale.pl    =          2; % path-loss factor
%     channel: antenna gain pattern (isotropic)
stdsettings.channel{1,1}.directivity.txant  =     '';
stdsettings.channel{1,1}.directivity.txrot  = [0, 0];
stdsettings.channel{1,1}.directivity.rxant  =     '';
stdsettings.channel{1,1}.directivity.rxrot  = [0, 0];
stdsettings.channel{1,1}.directivity.dir_tx = [0, 0];
%     channel: smallscale (OFF)
stdsettings.channel{1,1}.smallscale.on        = false;
stdsettings.channel{1,1}.smallscale.det       =   NaN;
stdsettings.channel{1,1}.smallscale.seed      =   NaN;
stdsettings.channel{1,1}.smallscale.ensembles =   NaN;
stdsettings.channel{1,1}.smallscale.maxiter   =   NaN;
stdsettings.channel{1,1}.smallscale.maxrays   =   NaN;
stdsettings.channel{1,1}.smallscale.k         =   NaN;
stdsettings.channel{1,1}.smallscale.trms      =   NaN;
stdsettings.channel{1,1}.smallscale.bw        =   NaN;
stdsettings.channel{1,1}.smallscale.fres      =   NaN;
stdsettings.channel{1,1}.smallscale.eps_k     =   NaN;
stdsettings.channel{1,1}.smallscale.eps_trms  =   NaN;
stdsettings.channel{1,1}.smallscale.fs        =   NaN;
stdsettings.channel{1,1}.noise.type =     'off';
stdsettings.channel{1,1}.noise.n0   = [NaN,NaN];
stdsettings.channel{1,1}.noise.frxs =       NaN;
stdsettings.channel{1,1}.noise.fs   =       NaN;       
%     mfcw_addseccarriers
stdsettings.mfcw_addseccarriers.range_s = []; % areas of signal to add secondary carriers to; leave empty to select entire signal
%     mfcw_compsel
stdsettings.mfcw_compsel.nharm  =   20; % number of harmonics to consider (for check only) 
stdsettings.mfcw_compsel.iirord =    4; % selectionbandpass filter order
stdsettings.mfcw_compsel.att    =  100; % stopband attenuation of bandpasses in dB
%     mfcw_calcdist
stdsettings.mfcw_calcdist.c_ord = []; % default carrier permutation

% random setups [min, max]
%     general
randsettings.fs   = [  2e9,   3e9]; % Hz sampling frequency
randsettings.frs  = [    3,     7]; % times max(fi)+fm reader sampling frequency
randsettings.fc   = [250e6, 500e6]; % Hz center carrier frequency
randsettings.fm   = [ 25e3,  50e3]; % Hz modulation frequency
randsettings.dist = [    1,     9]; % m distance (make sure this is within first period of fi for c below)
randsettings.c    = [  5e8,   3e9]; % m/s speed of light
randsettings.drho = [  0.1,   0.8]; % modulation depth (reflection coefficient max-min)
%     secondary carriers
randsettings.nc     = [    1,     5]; % number of secondary carriers to add (<= .ni)
randsettings.ni     = [    1,     5]; % length of secondary carrier definitions (>= .nc)
randsettings.fi     = [  1e6,  10e6]; % Hz secondary carrier frequency offsets
randsettings.vari   = [ 1e-3,   1e1]; % W secondary carrier variance


% *******************************************************************************************************
% check setup and complete/check partitioning 

if max(randsettings.dist) > min(randsettings.c) / (4*max(randsettings.fi)) 
   error('Distance/frequency/speed-of-light setup will lead to ambiguities.')
end

% partitioning
if length(partitioning.names) ~= length(partitioning.runs)
   error('Length of entries of struct partitioning do not match. Check partitioning.');
end
%     indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% create settings

for i = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)

      % start from standard setup
      settings{j} = stdsettings;

      % general settings
      settings{j}.fs   = randsettings.fs  (1) + rand*(randsettings.fs  (2) - randsettings.fs  (1));
      settings{j}.fc   = randsettings.fc  (1) + rand*(randsettings.fc  (2) - randsettings.fc  (1));
      settings{j}.fm   = randsettings.fm  (1) + rand*(randsettings.fm  (2) - randsettings.fm  (1));
      settings{j}.dist = randsettings.dist(1) + rand*(randsettings.dist(2) - randsettings.dist(1));
      settings{j}.c    = randsettings.c   (1) + rand*(randsettings.c   (2) - randsettings.c   (1));
      settings{j}.drho = randsettings.drho(1) + rand*(randsettings.drho(2) - randsettings.drho(1)); 
      %     round distance to sampling interval
      settings{j}.dist = round(settings{j}.dist/settings{j}.c*settings{j}.fs)*settings{j}.c/settings{j}.fs;
 
      % secondary carriers
      ni = intrand(randsettings.ni);
      settings{j}.nc = Inf;
      while settings{j}.nc > ni
         settings{j}.nc = intrand(randsettings.nc);
      end
      settings{j}.fi = cumsum( randsettings.fi(1) + rand(ni,1)*(randsettings.fi(2) - randsettings.fi(1)) );
      
      % needed vector length
      settings{j}.len_s = settings{j}.fs / settings{j}.fm * settings{j}.fm_per;
      
      % complete settings
      settings{j}.frs  = (max(settings{j}.fi) + settings{j}.fm) *...
         (randsettings.frs (1) + rand*(randsettings.frs (2) - randsettings.frs (1)));
      settings{j}.mfcw_addseccarriers.vari   = ...
         randsettings.vari(1) + rand(ni,1)*(randsettings.vari(2) - randsettings.vari(1));
      settings{j}.mfcw_addseccarriers.gammai = zeros(size(settings{j}.mfcw_addseccarriers.vari));
      
      % distribute settings
      settings{j}.channel{1,1}.largescale.fs   = settings{j}.fs;
      settings{j}.channel{1,1}.largescale.f0   = settings{j}.fc;
      settings{j}.channel{1,1}.largescale.c    = settings{j}.c;
      settings{j}.channel{1,1}.largescale.dist = settings{j}.dist;
      settings{j}.mfcw_addseccarriers.fs = settings{j}.fs;
      settings{j}.mfcw_addseccarriers.f0 = settings{j}.fc;
      settings{j}.mfcw_addseccarriers.nc     = settings{j}.nc;
      settings{j}.mfcw_addseccarriers.fi     = settings{j}.fi;
      settings{j}.mfcw_compsel.nc  = settings{j}.nc;
      settings{j}.mfcw_compsel.fi  = settings{j}.fi;
      settings{j}.mfcw_compsel.lf  = settings{j}.fm;
      settings{j}.mfcw_compsel.bw  = settings{j}.fm / 2;
      settings{j}.mfcw_compsel.frs = settings{j}.frs;
      settings{j}.mfcw_calcdist.nc = settings{j}.nc;
      settings{j}.mfcw_calcdist.fi = settings{j}.fi; 
      settings{j}.mfcw_calcdist.c  = settings{j}.c;
      
   end
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
characteristic = 'settings and results for selftest: mfcw_addseccarriers/compsel/calcdist.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'partitioning');
disp(sprintf('\nSettings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));

