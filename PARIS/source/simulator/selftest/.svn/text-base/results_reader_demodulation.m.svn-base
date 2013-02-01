% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_reader_demodulation function
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
filename  = 'results_reader_demodulation';

% nonrandom test setup
%     symbol constellations (data itself is random)
i = complex(0,1);
testsettings.con = [1+i, -1+i, -1-i, 1-i]/sqrt(2); % unit energy!
%     length of test vector if no data is modulated
testsettings.veclen = 2e6; % samples @ fs (should be significantly longer @fs than impres of demod filter)
%     oscillator (mixer) settings: clean carrier (otherwise noise problems are possible)
testsettings.fstddev =   0; % Hz frequency standard deviation (Gaussian) ... phase noise
testsettings.astddev =   0; % amplitude standard deviation (center = 1; Gaussian)
testsettings.snr     = Inf; % SNR in dBc (to carrier)

% random test setup bounds (uniformly distributed)
%     frequencies
randsettings.fs  = [1e9, 2e9]; % Hz sampling frequency
randsettings.f0  = [1e8, 2e8]; % Hz carrier frequency (clean cosine)
randsettings.frs = [5e7, 8e7]; % Hz reader sampling frequency
randsettings.fm  = [1e5, 3e5]; % Hz modulation frequency (rectangular modulation); rounded to frs
%     alti-aliasing-filter (IIR)
randsettings.iirord = [ 3,  6]; % filter order
randsettings.att    = [30, 50]; % dB stopband attenuation
randsettings.cutoff = [18, 24]; % times modulation frequency stopband edge frequency (integer!)
%     quantization
randsettings.q = [1, 16]; % bits quantization
%     lenghts
randsettings.len  = [0.9, 1.1]; % length factor (for truncation/zero-padding)
randsettings.symb = [ 31,  99]; % number of data symbols
%     other
randsettings.varc = [1e-1, 1e-1]; % carrier variance
randsettings.varn = [1e-3, 1e-2]; % noise variance

% partitioning
partitioning.names = {'noise', 'sinusoid', 'data', 'quantization'};
partitioning.runs  = [20, 25, 100, 50];


% *******************************************************************************************************
% check setup

% frequencies
%     OSR for f0
if randsettings.fs(1) < 5 * randsettings.f0(2)
   error('Low oversampling rate: f0 vs fs');
end
%     OSR for fcut (includes fm vs frs)
if randsettings.frs(1)/2 < randsettings.fm(2)*randsettings.cutoff(2) 
   error('Low oversampling rate: fcut vs frs');
end
%     excessive RAM usage
if randsettings.fs(2) / randsettings.fm(1) * randsettings.symb(2) > 1e7
   disp('Warning: High RAM usage');
end
if randsettings.fs(2) / randsettings.fm(1) > 1e8
   error('RAM usage too high');
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
           
      % sampling frequencies
      settings{j}.fs  = randsettings.fs (1) + rand*(randsettings.fs (2) - randsettings.fs (1)); % global
      settings{j}.frs = randsettings.frs(1) + rand*(randsettings.frs(2) - randsettings.frs(1)); % reader
            
      % "standard" definitions (might be overwritten below!)
      %     carrier and modulation frequency
      settings{j}.f0  = randsettings.f0 (1) + rand*(randsettings.f0 (2) - randsettings.f0 (1));
      settings{j}.fm  = randsettings.fm (1) + rand*(randsettings.fm (2) - randsettings.fm (1));
      %     quantization
      settings{j}.q = NaN; % no quantization
            
      % test specific mode 
      switch lower(partitioning.names{i})
         case 'noise'
            settings{j}.mode = struct('carrier',0, 'data',0, 'noise',1);
         case 'sinusoid'
            settings{j}.mode = struct('carrier',1, 'data',0, 'noise',0);
         case 'data'
            settings{j}.mode = struct('carrier',1, 'data',1, 'noise',0);
         case 'quantization'
            settings{j}.mode = struct('carrier',1, 'data',1, 'noise',1);
            settings{j}.q  = intrand(randsettings.q);
         otherwise
            error('settings generation: check partitioning');
      end
      
      % input signal settings
      %     carrier?
      if settings{j}.mode.carrier
         settings{j}.varc = randsettings.varc(1) + rand*(randsettings.varc(2) - randsettings.varc(1));
      else
         settings{j}.varc = 0;
      end
      %     data?
      if settings{j}.mode.data
         settings{j}.symb = intrand(randsettings.symb); % number of symbols
         settings{j}.ind_data = intrand(1,settings{j}.symb, [1,length(testsettings.con)]);
         settings{j}.data = testsettings.con(settings{j}.ind_data);
      else
         settings{j}.symb = 0;
         settings{j}.ind_data = [];
         settings{j}.data = [];
      end
      %     noise?
      if settings{j}.mode.noise
         settings{j}.varn = randsettings.varn(1) + rand*(randsettings.varn(2) - randsettings.varn(1));
      else
         settings{j}.varn = 0;
      end
                    
      % setup of anti-aliasing filter (IIR)
      settings{j}.iirord = intrand(randsettings.iirord);
      settings{j}.att  = randsettings.att(1) + rand*(randsettings.att(2) - randsettings.att(1));
      settings{j}.fcut = settings{j}.fm * (randsettings.cutoff(1) + ...
                         rand*(randsettings.cutoff(2) - randsettings.cutoff(1)));  
      
      % set lengths ("automatic" for data)
      if settings{j}.mode.data
         % input vector length (... better too long than too short)
         settings{j}.veclen_s = settings{j}.symb * ceil(settings{j}.fs / settings{j}.fm);
         % output vector length
         settings{j}.len_rs = NaN; % automatic
               
      else
         % input vector length
         settings{j}.veclen_s = testsettings.veclen;
         % output vector length (will lead to critical warnings)
         settings{j}.len_rs = round(settings{j}.veclen_s * settings{j}.frs / settings{j}.fs *... % approx length @ frs
            (randsettings.len(1)+rand*(randsettings.len(2)-randsettings.len(1))));
      end
      
      % prepare mixer (oscillator) settings
      oscsettings{j}.fcenter = settings{j}.f0; % Hz center frequency
      oscsettings{j}.fstddev = testsettings.fstddev; % Hz frequency standard deviation (Gaussian) ... phase noise
      oscsettings{j}.astddev = testsettings.astddev; % amplitude standard deviation (center = 1; Gaussian)
      oscsettings{j}.snr     = testsettings.snr; % SNR in dBc (to carrier)
      oscsettings{j}.length  = settings{j}.veclen_s/settings{j}.fs; % s length of carrier wave signal
      oscsettings{j}.mode    = 'udef'; % will be set by reader_demodulation
      oscsettings{j}.fs      = settings{j}.fs; % Hz sampling frequency
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
characteristic = 'settings and results for selftest: reader_demodulation.m'; %#ok<NASGU>
createdon      = datestr(now, 0); %#ok<NASGU>
createdby      = sprintf('%s.m, rev.: %s', mfilename, version); %#ok<NASGU>
usedmatfiles   = ''; %#ok<NASGU>

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'oscsettings', 'partitioning');
disp(sprintf('Settings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));
