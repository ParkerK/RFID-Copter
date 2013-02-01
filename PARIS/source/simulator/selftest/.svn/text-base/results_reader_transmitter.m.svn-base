% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_reader_transmitter function
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
filename  = 'results_reader_transmitter';

% nonrandom test setup
%     periodogram estimator setup
testsettings.ns = 64; % samples between stopband edge frequencies => resolution
%     minimum number of averaged spectra (NFFT*?-NFFT/2 <= vector length)
testsettings.minspec = 6;

% random test setup bounds (uniformly distributed)
%     basic settings
randsettings.fs  = [ 1e8, 5e8]; % Hz sampling frequency
randsettings.len = [ 1e5, 1e6]; % samples @ fs vector length
randsettings.var = [1e-3, 1e1]; % input signal variance
randsettings.ptx = [1e-2, 1e2]; % output (transmit) variance
%     bandpass (if enabled)
randsettings.bp_ord = [ 4,  6]; % filter order
randsettings.bp_att = [20, 40]; % dB stopband attenuation
randsettings.bp_fc  = [5e6, 1e7]; % Hz (logarithmic) center frequency
randsettings.bp_bw  = [1e5, 1e6]; % Hz bandwidth
%     nonlinear poweramp (if enabled)
randsettings.nl_lim = [0.1, 1]; % symmetrical signal amplitude limits 0<=x<=1 (1: maximum limits of pwramp)
randsettings.nl_char = {'readerchar_transmitter'}; % filenames for nonlinear power amplifier characteristic

% partitioning
partitioning.names = {'linear amp', 'linear amp + bandpass', 'pwramp', 'pwramp + bandpass'};
partitioning.runs  = [10, 25, 25, 25];

% tolerances for reverse checks
internalsettings.tol_fc = 1; % Hz tolerance for filter center frequency
internalsettings.tol_bw = 1; % Hz tolerance for filter bandwidth


% *******************************************************************************************************
% check setup (simple checks => false positives/negatives possible)

% make sure stopband attuenuation can be reached (assuming 20 dB/dec for n-th order filter)
%     f <<
if 20*randsettings.bp_ord(1)*log10(randsettings.bp_fc(1)-randsettings.bp_bw(2)) < randsettings.bp_att(2) % @ 1 Hz
   error('Stopband attuenuation cannot be guaranteed for small frequencies. Check setup.');
end
%     f >>
if 20*randsettings.bp_ord(1)*(log10(randsettings.fs(1)/2)-log10(randsettings.bp_fc(2)+randsettings.bp_bw(2))) < randsettings.bp_att(2)
   error('Stopband attuenuation cannot be guaranteed for high frequencies. Check setup.');
end

% make sure 0 dB passband frequency can be reached (assuming 20 dB/dec for n-th order filter)
if 20*randsettings.bp_ord(1)*log10(randsettings.bp_bw(1)) < randsettings.bp_att(1)
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
      settings{j}.fs      = randsettings.fs (1) + rand*(randsettings.fs (2) - randsettings.fs (1));
      settings{j}.var     = randsettings.var(1) + rand*(randsettings.var(2) - randsettings.var(1));
      settings{j}.len     = intrand(randsettings.len);

      % bandpass settings
      %     filter order / attenuation
      settings{j}.bp_ord = intrand(randsettings.bp_ord);
      settings{j}.bp_att  = randsettings.bp_att(1) + rand*(randsettings.bp_att(2) - randsettings.bp_att(1));
      %     cutoff via (logarithmic) center frequency and bandwidth
      fc = randsettings.bp_fc(1) + rand*(randsettings.bp_fc(2) - randsettings.bp_fc(1));
      bw = randsettings.bp_bw(1) + rand*(randsettings.bp_bw(2) - randsettings.bp_bw(1));
      settings{j}.bp_fcut = [(sqrt(bw^2+4*fc^2)-bw)/2, 2*fc^2/(sqrt(bw^2+4*fc^2)-bw)]; % [lower, upper]
      %     reverse check for cutoff
      if abs(settings{j}.bp_fcut(2)-settings{j}.bp_fcut(1)-bw) >= internalsettings.tol_bw ||... % check bandwidth
         abs(settings{j}.bp_fcut(2)*settings{j}.bp_fcut(1)-fc^2) >= internalsettings.tol_fc % check center frequency
         error('filter reverse check failed: check calculation');
      end
      
      % amplifier settings
      settings{j}.ptx     = randsettings.ptx(1)    + rand*(randsettings.ptx(2)    - randsettings.ptx(1)   );
      settings{j}.nl_lim  = randsettings.nl_lim(1) + rand*(randsettings.nl_lim(2) - randsettings.nl_lim(1));
      settings{j}.nl_char = randsettings.nl_char{ intrand([1, length(randsettings.nl_char)]) };
      
      % calculate needed window length for spectrum estimator to meet given resolution
      settings{j}.nfft = round(testsettings.ns * settings{j}.fs/2 / diff(settings{j}.bp_fcut));
      settings{j}.nfft = 2^ceil(log2(settings{j}.nfft)); % power of 2 for fft
      
      % check FFT length vs. vector length
      if settings{j}.nfft * testsettings.minspec - settings{j}.nfft/2 > settings{j}.len
         disp(sprintf('   Warning: Insufficient vector length. Increasing from %i to %i samples.', ...
            settings{j}.len, settings{j}.nfft * testsettings.minspec - settings{j}.nfft/2));
         settings{j}.len = settings{j}.nfft * testsettings.minspec - settings{j}.nfft/2;
      end

      % test specific settings 
      switch lower(partitioning.names{i})
         case 'linear amp'
            settings{j}.en_nl = false;
            settings{j}.en_bp = false;
         case 'linear amp + bandpass'
            settings{j}.en_nl = false;
            settings{j}.en_bp = true;
         case 'pwramp'
            settings{j}.en_nl = true;
            settings{j}.en_bp = false;
         case 'pwramp + bandpass'
            settings{j}.en_nl = true;
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
characteristic = 'settings and results for selftest: reader_transmitter.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
