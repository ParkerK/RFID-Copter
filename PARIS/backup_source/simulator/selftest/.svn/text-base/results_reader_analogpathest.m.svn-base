% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_reader_analogpathest function
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
filename  = 'results_reader_analogpathest';

% nonrandom test setup
%     general
stdsettings.ens        = 1e5; % number of ensembles for random estimation error tests
stdsettings.fieldnames = {'transmitter'; 'receiver'; 'demodulation'; 'all'}; % field names of returned struct

%     reader_transmitter
stdsettings.reader_transmitter.en_nl   = false; % pwramp: enable nonlinear power amplifier
stdsettings.reader_transmitter.nl_char =    ''; % pwramp: filename for nonlinear power amplifier characteristic
stdsettings.reader_transmitter.nl_lim  =   NaN; % pwramp: symmetrical input signal amplitude limits 0<=x<=1 (1: max limits of pwramp)
%     reader_demodulation
stdsettings.reader_demodulation.frs    = NaN; % Hz reader sampling frequency
stdsettings.reader_demodulation.q      = NaN; % bits quantization (set to NaN to deactivate)
stdsettings.reader_demodulation.len_rs = NaN; % length of output signals in samples @ frs (NaN to deactivate)
%     reader_oscillator (just a dummy for reader_demodulation)
stdsettings.reader_oscillator.fcenter = NaN;
stdsettings.reader_oscillator.fstddev = NaN;
stdsettings.reader_oscillator.astddev = NaN;
stdsettings.reader_oscillator.snr     = NaN;
stdsettings.reader_oscillator.length  = NaN;
stdsettings.reader_oscillator.fs      = NaN;
stdsettings.reader_oscillator.mode    =  '';

% random test setup bounds (uniformly distributed) [min, max]
%     general
randsettings.fs       = [ 5e9, 27e9]; % Hz sampling frequency
randsettings.fc       = [ 5e7,  2e9]; % Hz carrier center frequency
randsettings.freq_len = [   1,   10]; % length of frequency vector (test will use [-freq, freq])
randsettings.freq     = [-1e7,  1e7]; % Hz baseband frequency bounds (make sure this can be outside passbands)
%     reader_transmitter
randsettings.reader_transmitter.p_enbp =         0.5; % probability for enabled bandpass
randsettings.reader_transmitter.ptx    = [1e-2, 1e2]; % output (transmit) variance
randsettings.reader_transmitter.bp_ord = [   4,   6]; % filter order
randsettings.reader_transmitter.bp_att = [  20,  40]; % dB stopband attenuation
randsettings.reader_transmitter.bp_fc  = [ 5e6, 1e7]; % Hz (logarithmic) center frequency
randsettings.reader_transmitter.bp_bw  = [ 1e5, 1e6]; % Hz bandwidth
%     reader_demodulation
randsettings.reader_demodulation.p_enbp =       0.5; % probability for enabled bandpass
randsettings.reader_demodulation.iirord = [  3,   6]; % filter order
randsettings.reader_demodulation.att    = [ 30,  50]; % dB stopband attenuation
randsettings.reader_demodulation.cutoff = [0.1, 0.4]; % times sampling frequency
%     reader_receiver
randsettings.reader_receiver.p_enbp =        0.5; % probability for enabled bandpass
randsettings.reader_receiver.iirord = [  4,   6]; % filter order
randsettings.reader_receiver.att    = [ 20,  40]; % dB stopband attenuation
randsettings.reader_receiver.fc     = [5e6, 1e7]; % Hz (logarithmic) center frequency
randsettings.reader_receiver.bw     = [1e5, 1e6]; % Hz bandwidth
%     reader_analogpathest
randsettings.reader_analogpathest.err_syst = [1/sqrt(2),      pi]; % systemativ gain
randsettings.reader_analogpathest.err_rand = [     1/pi, sqrt(2)]; % random gain

% partitioning
partitioning.names = {'ideal', 'systematic errors', 'random errors', 'systematic and random errors'};
partitioning.runs  = [250, 250, 10, 10];


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
      
      % copy from stdsettings
      settings{j} = stdsettings;
      
      % general settings
      %     frequency settings
      settings{j}.fs   =       randsettings.fs  (1) + rand*(randsettings.fs  (2) - randsettings.fs  (1));
      settings{j}.fc   =       randsettings.fc  (1) + rand*(randsettings.fc  (2) - randsettings.fc  (1));
      %     distribute
      settings{j}.reader_transmitter.fs   = settings{j}.fs;
      settings{j}.reader_transmitter.fc   = settings{j}.fc;
      settings{j}.reader_demodulation.fs  = settings{j}.fs;
      settings{j}.reader_receiver.fs      = settings{j}.fs;
      settings{j}.reader_analogpathest.fs = settings{j}.fs;
      settings{j}.reader_analogpathest.f0 = settings{j}.fc;
      %     checked frequencies
      settings{j}.freq_len = round(randsettings.freq_len(1) +...
         rand*(randsettings.freq_len(2) - randsettings.freq_len(1)));
      settings{j}.freq = round(randsettings.freq(1) +...
         rand(settings{j}.freq_len,1)*(randsettings.freq(2) - randsettings.freq(1)));
      %     only one ensemble for nonrandom tests
      if ~strcmpi(partitioning.names(i), 'random errors') &&...
            ~strcmpi(partitioning.names(i), 'systematic and random errors')
         settings{j}.ens = 1;
      end
      
      % settings for reader_transmitter
      %     general settings
      settings{j}.reader_transmitter.en_bp = round(randsettings.reader_transmitter.p_enbp * rand);
      %     bandpass settings
      settings{j}.reader_transmitter.bp_ord = intrand(randsettings.reader_transmitter.bp_ord);
      settings{j}.reader_transmitter.bp_att = randsettings.reader_transmitter.bp_att(1) +...
         rand*(randsettings.reader_transmitter.bp_att(2) - randsettings.reader_transmitter.bp_att(1));
      fc = randsettings.reader_transmitter.bp_fc(1) + ...
         rand*(randsettings.reader_transmitter.bp_fc(2) - randsettings.reader_transmitter.bp_fc(1));
      bw = randsettings.reader_transmitter.bp_bw(1) + ...
         rand*(randsettings.reader_transmitter.bp_bw(2) - randsettings.reader_transmitter.bp_bw(1));
      settings{j}.reader_transmitter.bp_fcut = [(sqrt(bw^2+4*fc^2)-bw)/2, 2*fc^2/(sqrt(bw^2+4*fc^2)-bw)]; % [lower, upper]
      %     (linear) amplifier settings
      settings{j}.reader_transmitter.ptx = randsettings.reader_transmitter.ptx(1) +...
            rand*(randsettings.reader_transmitter.ptx(2) - randsettings.reader_transmitter.ptx(1));

      % settings for reader_demodulation
      settings{j}.reader_demodulation.iirord = intrand(randsettings.reader_demodulation.iirord);
      settings{j}.reader_demodulation.att  = randsettings.reader_demodulation.att(1) +...
         rand*(randsettings.reader_demodulation.att(2) - randsettings.reader_demodulation.att(1));
      settings{j}.reader_demodulation.fcut = settings{j}.fs * (randsettings.reader_demodulation.cutoff(1) + ...
         rand*(randsettings.reader_demodulation.cutoff(2) - randsettings.reader_demodulation.cutoff(1)));
      
      % settings for reader_receiver
      %     general settings
      settings{j}.reader_receiver.en_bp = round(randsettings.reader_receiver.p_enbp * rand);
      %     bandpass settings
      settings{j}.reader_receiver.iirord = intrand(randsettings.reader_receiver.iirord);
      settings{j}.reader_receiver.att = randsettings.reader_receiver.att(1) +...
         rand*(randsettings.reader_receiver.att(2) - randsettings.reader_receiver.att(1));
      fc = randsettings.reader_receiver.fc(1) +...
         rand*(randsettings.reader_receiver.fc(2) - randsettings.reader_receiver.fc(1));
      bw = randsettings.reader_receiver.bw(1) +...
         rand*(randsettings.reader_receiver.bw(2) - randsettings.reader_receiver.bw(1));
      settings{j}.reader_receiver.fcut = [(sqrt(bw^2+4*fc^2)-bw)/2, 2*fc^2/(sqrt(bw^2+4*fc^2)-bw)]; % [lower, upper]
      
      % settings.for reader_analogpathest
      %     systematic gain
      if strcmpi(partitioning.names(i), 'systematic errors') ||...
            strcmpi(partitioning.names(i), 'systematic and random errors')
         settings{j}.reader_analogpathest.err(1) = randsettings.reader_analogpathest.err_syst(1) +...
            rand*(randsettings.reader_analogpathest.err_syst(2) - randsettings.reader_analogpathest.err_syst(1));
      else
         settings{j}.reader_analogpathest.err(1) = 1; % no systematic error => systematic gain = 1
      end
      %     random gain
      if strcmpi(partitioning.names(i), 'random errors') ||...
            strcmpi(partitioning.names(i), 'systematic and random errors')
         settings{j}.reader_analogpathest.err(2) = randsettings.reader_analogpathest.err_rand(1) +...
            rand*(randsettings.reader_analogpathest.err_rand(2) - randsettings.reader_analogpathest.err_rand(1));
      else
         settings{j}.reader_analogpathest.err(2) = 0; % no random error => stddev of gain = 0
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
characteristic = 'settings and results for selftest: reader_analogpathest.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
