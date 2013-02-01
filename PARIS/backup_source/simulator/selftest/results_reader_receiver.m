% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_reader_receiver function
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
globalinit('silent');

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
filename  = 'results_reader_receiver';

% nonrandom test setup
%     periodogram estimator setup
testsettings.ns = 64; % samples between stopband edge frequencies => resolution
%     minimum number of averaged spectra (NFFT*?-NFFT/2 <= vector length)
testsettings.minspec = 6;

% random test setup bounds (uniformly distributed)
%     basic settings
randsettings.fs  = [ 1e8, 5e8]; % Hz sampling frequency
randsettings.len = [ 1e5, 1e6]; % samples @ fs vector length
randsettings.var = [1e-2, 1e2]; % input signal variance
%     bandpass (if enabled)
randsettings.iirord = [ 4,  6]; % filter order
randsettings.att    = [20, 40]; % dB stopband attenuation
randsettings.fc     = [5e6, 1e7]; % Hz (logarithmic) center frequency
randsettings.bw     = [1e5, 1e6]; % Hz bandwidth

% expected values
%     center attenuation
results.att_fc = -6; % dB (power splitter has factor 1/2)

% partitioning
partitioning.names = {'power splitter', 'bandpass + power splitter'};
partitioning.runs  = [ 5, 25];

% tolerances for reverse checks
internalsettings.tol_fc = 1; % Hz tolerance for filter center frequency
internalsettings.tol_bw = 1; % Hz tolerance for filter bandwidth


% *******************************************************************************************************
% check setup (simple checks => false positives/negatives possible)

% make sure stopband attuenuation can be reached (assuming 20 dB/dec for n-th order filter)
%     f <<
if 20*randsettings.iirord(1)*log10(randsettings.fc(1)-randsettings.bw(2)) < randsettings.att(2) % @ 1 Hz
   error('Stopband attuenuation cannot be guaranteed for small frequencies. Check setup.');
end
%     f >>
if 20*randsettings.iirord(1)*(log10(randsettings.fs(1)/2)-log10(randsettings.fc(2)+randsettings.bw(2))) < randsettings.att(2)
   error('Stopband attuenuation cannot be guaranteed for high frequencies. Check setup.');
end

% make sure 0 dB passband frequency can be reached (assuming 20 dB/dec for n-th order filter)
if 20*randsettings.iirord(1)*log10(randsettings.bw(1)) < randsettings.att(1)
   error('Passband attenuation of 0dB cannot be guaranteed. Check setup.');
end


% *******************************************************************************************************
% complete/check partitioning

if length(partitioning.names) ~= length(partitioning.runs)
   criterr('Length of entries of struct partitioning do not match. Check partitioning.');
end

% indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% create settings
disp(sprintf('Creating results and settings.'));

for i = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(i) : partitioning.indices(i+1)
      
      % basic settings
      settings{j}.fs  = randsettings.fs (1) + rand*(randsettings.fs (2) - randsettings.fs (1));
      settings{j}.var = randsettings.var(1) + rand*(randsettings.var(2) - randsettings.var(1));
      settings{j}.len = intrand(randsettings.len);
      
      % bandpass settings
      %     filter order / attenuation
      settings{j}.iirord = intrand(randsettings.iirord);
      settings{j}.att    = randsettings.att(1) + rand*(randsettings.att(2) - randsettings.att(1));
      %     cutoff via (logarithmic) center frequency and bandwidth
      fc = randsettings.fc(1) + rand*(randsettings.fc(2) - randsettings.fc(1));
      bw = randsettings.bw(1) + rand*(randsettings.bw(2) - randsettings.bw(1));
      settings{j}.fcut = [(sqrt(bw^2+4*fc^2)-bw)/2, 2*fc^2/(sqrt(bw^2+4*fc^2)-bw)]; % [lower, upper]
      %     reverse check for cutoff
      if abs(settings{j}.fcut(2)-settings{j}.fcut(1)-bw) >= internalsettings.tol_bw ||... % check bandwidth
         abs(settings{j}.fcut(2)*settings{j}.fcut(1)-fc^2) >= internalsettings.tol_fc % check center frequency
         error('filter reverse check failed: check calculation');
      end
      
      % calculate needed window length for spectrum estimator to meet given resolution
      settings{j}.nfft = round(testsettings.ns * settings{j}.fs/2 / diff(settings{j}.fcut));
      settings{j}.nfft = 2^ceil(log2(settings{j}.nfft)); % power of 2 for fft
      
      % check FFT length vs. vector length
      if settings{j}.nfft * testsettings.minspec - settings{j}.nfft/2 > settings{j}.len
         disp(sprintf('   Warning: Insufficient vector length. Increasing from %i to %i samples.', ...
            settings{j}.len, settings{j}.nfft * testsettings.minspec - settings{j}.nfft/2));
         settings{j}.len = settings{j}.nfft * testsettings.minspec - settings{j}.nfft/2;
      end
                  
      % test specific settings 
      switch i
         case 1 % only power splitter
            settings{j}.en_bp = false;
         case 2 % bandpass plus power splitter
            settings{j}.en_bp = true;
         otherwise
            error('settings generation: check partitioning');
      end
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
characteristic = 'settings and results for selftest: reader_receiver.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'results', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
