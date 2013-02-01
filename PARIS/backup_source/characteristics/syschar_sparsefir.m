% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate data for: system - sparse FIR implementation
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


% ******************************************************************************************************
% initialization
clear; close all; clc;
version = 'beta 3.0';

% one last chance to reconsider
disp(sprintf('WARNING: This will overwrite existing files and might take a few days!\n'));
reply = input('Are you REALLY sure y/n [n]?  ', 's');
if ~strcmpi(reply, 'y')
   return
end
disp(' ');

% paths, etc
%     path to globalinit.m
path = fullfile(fileparts(mfilename('fullpath')), '..');
%     Matlab does not resolve symbolic links; if we're running on Unix, let the system resolve them
if isunix
   [dummy, path] = system(sprintf('readlink -f "%s"', path));
end
%     add this path and call globalinit script
restoredefaultpath; addpath(path); clear('dummy', 'path');
globalsettings_copy = globalinit('exceptions');

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% check if hostname is ok
if isempty(globalsettings_copy.misc.hostname)
   disp('WARNING: hostname is empty; filename will end with _');
end


% ******************************************************************************************************
% settings

% number of data points
% ... performance of interp1 should be roughly constant up to a few hundred points
settings.n_i =  25; % # input signal lengths
settings.n_f =  20; % # filter length
settings.n_z =  25; % maximum # of test points between sparsity 0 and 1 (settings.eta_f depends on this value!)

% bounds (powers of 10)
settings.exp_i = [       2, 7]; % [start, stop] exponent for input signal lengths (log10)
settings.exp_f = [log10(2), 3]; % [start, stop] exponent for filter lengths (log10) d=10m @ fs=30GHz: 1000 taps

% interpolation to get tipping point and final filtering
settings.n_int       = 1e2; % # of pts for interpolation sparsity 0...1 to get tipping point
settings.medfilt_len =   3; % smoothing filter for zero-crossings and final characteristic

% average seconds time per sample of input vector and filter vector for ETA estimate (n*logn) 
% ... depends on host (only valid for high signal lengths) and settings.n_z
% ... only valid for long execution times; exp_i=[2,7], exp_f=[log10(2), 3], n_i/f/z=25/20/25
switch lower(globalsettings_copy.misc.hostname)
   case 'wishmaster'
      settings.eta_f = 1.8e-8;
   case 'apsalar'
      settings.eta_f = 2.9e-8;
   case {'deralte', 'isis', 'aphrodite'}
      settings.eta_f = 3.1e-8;
   case 'edgewalker'
      settings.eta_f = 5.0e-8; % 5.25 blas
   case 'hammer01'
      settings.eta_f = 5.0e-8;
   case {'hammer02', 'hammer03'}
      settings.eta_f = 5.8e-8; % 2 nblas/blas: 5.84/5.75, 3 nblas: 5.66
   case 'werewolf'
      settings.eta_f = 6.8e-8; % ?
   case 'hammer05'
      disp('WARNING: Speed of hammer05 might have changed. ETA might be wrong.'); disp(' ');
      settings.eta_f = 2.1e-8; % 2.65 blas
   otherwise
      disp(sprintf('WARNING: Host "%s" not in speed estimate list. ETA might be off by factors!\n',...
         globalsettings_copy.misc.hostname));
      settings.eta_f = 2.0e-8;
end

% other
settings.numthreads    = 4; % maximum number of threads (will be saturated to number of CPUs available)
settings.renew_kticket = 0; % renew Kerberos ticket (this might take some time)?

% some checks
if min(settings.exp_i) < log10(2) || min(settings.exp_i) < log10(2)
   error('Lengths must be greater or equal two.');
end


% ******************************************************************************************************
% measure execution time of filter/loop implementation to get characteristic

% warn the user about the simple ETA calculation
disp('Note that the ETA ist only a very coarse estimate. If the end time starts to increase significantly');
disp('for larger loop values, expect a considerably longer simulation time.');

% multithreading
%     determine maximum amount of threads
maxNumCompThreads('automatic'); 
settings.maxthreads = maxNumCompThreads;
%     set to some value
maxNumCompThreads(min(settings.maxthreads, settings.numthreads));

% create length vectors
v_len_i = round( logspace(settings.exp_i(1), settings.exp_i(2), settings.n_i) ); % [samples] 
v_len_f = round( logspace(settings.exp_f(1), settings.exp_f(2), settings.n_f) ); % [samples]
%     make sure there are no identical entries
% v_len_i(diff(v_len_i)==0) = [];
% v_len_f(diff(v_len_f)==0) = [];
if any(diff(v_len_i)==0)
   error('Identical entries in input signal length vector (len_i). Check settings.');
end
if any(diff(v_len_f)==0)
   error('Identical entries in filter length vector (len_f). Check settings.');
end
%     "number of operations" for ETA
nop_eta = (v_len_i.*log10(v_len_i))' * (v_len_f.*log10(v_len_f));

% pre-initialize
results.sxing = nan(settings.n_i, settings.n_f);

% for length of input signal
for i = 1 : settings.n_i
   t_loopstart = clock;
   
   % renew Kerberos ticket and AFS authentication token if requested
   if settings.renew_kticket
      [status, msg] = renew_ticket();
      if status ~= 0
         disp(sprintf('   WARNING: Kerberos ticket and/or AFS authentication token renewal failed.'));
      end
      clear('status', 'msg');
   end
   
   % initialize
   len_i = v_len_i(i);
   input = randn(len_i,1);

   % for maximum filter length
   for f = 1 : settings.n_f       
      % initialize
      len_f = v_len_f(f);
      delay_full = [0 : 1 : len_f - 1]';
      gain_full  = 1 + rand(len_f, 1) + eps; % do not allow zeros
      
      % for different numbers of zeros
      v_z = round(linspace(1, len_f, min(len_f, settings.n_z)));
      if any(diff(v_z)==0)
         disp('   WARNING: Not strictly monotonic v_z found.');
      end
      v_z(diff(v_z)==0) = []; % make sure this is strictly monotonic (shouldn't happen anyway)
      % ... ATTENTION: zeros are inserted in a block instead of interleaved to avoid disambiguities
      for z = 1 : length(v_z)
         % record sparsity (number of zeros in filter)
         results.loop{i,f}.sparsity(z) = (v_z(z)-1) / len_f;
         
         % loop implementation
         %     initialize
         delay  = delay_full(v_z(z):end);
         gain   = gain_full (v_z(z):end);
         output1 = zeros(size(input));
         %     run
         tic;
         for k = 1 : length(delay)
            output1(1+delay(k):len_i) = output1(1+delay(k):len_i) + input(1:end-delay(k)) * gain(k);
         end
         results.loop{i,f}.t_loop(z) = toc;
         
         % filter implementation
         %     initialize filter (this is an additional step: take into account for time)
         tic;
         gain             = zeros(size(gain_full));
         gain(v_z(z):end) = gain_full(v_z(z):end);
         %     run
         output2 = filter(gain, 1, input);
         results.loop{i,f}.t_filter(z) = toc;
      end

      % error in case both methods created a different result
      if any( mean( (output1 - output2).^2 ) > eps)
         error('Outputs are not equal');
      end
      
      % find tipping point (decrease smoothing if zero-crossing is lost)
      for medfilt_len = settings.medfilt_len : -1 : 1
         %     interpolate to get better results
         s_int = linspace(0,1,settings.n_int);
         results.loop{i,f}.interp = interp1(results.loop{i,f}.sparsity,...
            medfilt1(results.loop{i,f}.t_filter ./ results.loop{i,f}.t_loop, medfilt_len), s_int, 'spline');
         %     find sparsity with identical exec. times
         ind = findzeros( results.loop{i,f}.interp - 1 );
         %     exit if there is a zero-crossing
         if ~isempty(ind); break; end
      end
      % record zero-crossing
      if ~isempty(ind)
         results.sxing(i,f) = s_int(ind(1)); % record only the first crossing
      else
         if all(results.loop{i,f}.interp < 1)
            results.sxing(i,f) = 1; % filter is always faster
         else
            results.sxing(i,f) = 0; % loop is always faster
         end
      end
   end
   
   % record time for manual ETA
   results.looptime(i) = etime(clock, t_loopstart);
   %     ETA based on "n*log(n)" assumption for filter and average loop_filter behavior (coarse)
   if i < settings.n_i
      results.eta(i) = settings.eta_f * sum(sum(nop_eta(i+1:end, :)));  
      eta_e = datevec(results.eta(i) / 86400); % ETA
      eta_d = feval('clock') + eta_e; % ETA: date/time
      disp(sprintf('* Loop %2i of %2i completed   -   ETA:  %2i days %2i hrs %2i min %2i sec (End: %s)',...
         i, settings.n_i, round(eta_e(3:end)), datestr(eta_d, 0)));
   else
      results.eta(i) = 0;
      eta_e = datevec(sum(results.looptime) / 86400); % ETA
      eta_d = feval('clock'); % overall time: date/time
      disp(sprintf('* Loop %2i of %2i completed   -   TIME: %2i days %2i hrs %2i min %2i sec (End: %s)',...
         i, settings.n_i, round(eta_e(3:end)), datestr(eta_d, 0)));
   end
end

% median filtering 
sxing = medfilt1(results.sxing', settings.medfilt_len)'; % along filter length
sxing = medfilt1(        sxing , settings.medfilt_len) ; % along input signal length

% average time for ETA estimate
disp(sprintf('*****\n ideal settings.eta_f = %e s', sum(results.looptime) / sum(sum(nop_eta)) ));

% cleanup
clear('input','output1','output2');


% ******************************************************************************************************
% save result

% prepare everything
%     header
matfilename    = sprintf('syschar_sparsefir_%s', globalsettings_copy.misc.hostname); %#ok<NASGU>
characteristic = sprintf('%s\n%s',...
                 'sparsity at performance crossing: sxing(length input, length filtercoeff)',...
                 'if filtercoeff sparsity < sxing (density > 1-sxing): filter() faster'); %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>
%     data (make sure all vectors are column-vectors for interp1q)
len_i       =       v_len_i(:);
len_i_log10 = log10(v_len_i(:));
len_f       =       v_len_f(:);
len_f_log10 = log10(v_len_f(:));
ind_i       = [1 : settings.n_i]';
ind_f       = [1 : settings.n_f]';
% save
save(fullfile('system',matfilename), 'matfilename','characteristic','createdby','usedmatfiles',...
   'len_i','len_i_log10','ind_i','len_f','len_f_log10','ind_f','sxing','settings');

% save entire workspace for debugging purposes and plots (it takes hours to run this script)
save(fullfile('workspaces', sprintf('%s_workspace', matfilename)));
