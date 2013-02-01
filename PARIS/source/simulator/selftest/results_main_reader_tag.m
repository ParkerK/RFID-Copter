% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% creates results for test_main_reader_tag function
% !!! EXISTING RESULTS WILL BE OVERWRITTEN - USE WITH CARE !!!
%
% This script is based upon sim_stdsettings_mfcw alpha 1.5, sim_mfcw_ranging alpha 1.2, and
% mfcw_statistics alpha 1.1.
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
% *******************************************************************************************************
% TEST SETUP AND RANDOM SETTINGS

% output directory and filename
directory = 'results';
filename  = 'results_main_reader_tag';
    
% partitioning
partitioning.names = {'data R->T', 'data T->R', 'MFCW'};
partitioning.runs  = [5, 5, 3];

% random setups [min, max]
%     general
randsettings.fs   = [  2e9,    3e9]; % Hz sampling frequency
randsettings.fc   = [860e6,  960e6]; % Hz center carrier frequency (keep close to UHF-RFID frequ. range)
%     reader
randsettings.frs  = [ 60e6,  80e6]; % Hz reader sampling frequency
randsettings.ptx  = [    3,     7]; % W transmit power (keep high-power)
randsettings.qy_s = {'s0', 's1', 's2', 's3'}; % query command: session
randsettings.qy_t =               {'a', 'b'}; % query command: target
randsettings.qy_q =                  [0, 16]; % query command: Q
%     tag
randsettings.fm   = [  24e3,   48e3]; % Hz modulation frequency (LF)
randsettings.pwrc = {'Vdda_V2A', 'Vdda_V2B', 'Vdda_V2C', 'Vdda_AVG'}; % power supply characteristics
randsettings.modc = {'tagchar_modulator_synth01', 'tagchar_modulator_synth02'}; % modulator characteristics
%     MFCW ranging
randsettings.dist = [  2.5,     5]; % m distance
randsettings.nc   = [    1,     3]; % number of secondary carriers (2FCW: nc=1); attention: nyquist with .frs
randsettings.fi   = [  2e6,   5e6]; % Hz secondary carrier frequency offsets; attention: nyquist with .frs
randsettings.vari = [ 1e-3,  1e-1]; % W secondary carrier variance (relative to var=1 for main carrier)



% *******************************************************************************************************
% *******************************************************************************************************
% STANDARD SETTINGS (based upon sim_stdsettings_mfcw alpha 1.5)

% *******************************************************************************************************
% general

stdsettings.c = 299792458; % m/s speed of light


% *******************************************************************************************************
% ranging

% multifrequency-CW-radar (mfcw)
%     secondary carrier generation
stdsettings.reader.ranging.mfcw_addseccarriers.nc   =         NaN; % number of secondary carriers
stdsettings.reader.ranging.mfcw_addseccarriers.fi   =          []; % Hz secondary carrier frequencies
stdsettings.reader.ranging.mfcw_addseccarriers.vari =          []; % secondary carrier variances
stdsettings.reader.ranging.mfcw_addseccarriers.f0   =         NaN; % Hz primary carrier frequency
stdsettings.reader.ranging.mfcw_addseccarriers.fs   =         NaN; % Hz sampling frequency
%     component selection
stdsettings.reader.ranging.mfcw_compsel.nc     =  NaN; % number of secondary carriers
stdsettings.reader.ranging.mfcw_compsel.fi     =   []; % Hz secondary carrier offset frequencies
stdsettings.reader.ranging.mfcw_compsel.lf     =  NaN; % Hz backscatter link frequency
stdsettings.reader.ranging.mfcw_compsel.bw     =  NaN; % +/- Hz stopband edge frequency
stdsettings.reader.ranging.mfcw_compsel.nharm  =   20; % number of harmonics to consider (for check only) 
stdsettings.reader.ranging.mfcw_compsel.iirord =    2; % selectionbandpass filter order
stdsettings.reader.ranging.mfcw_compsel.att    =  100; % stopband attenuation of bandpasses in dB
stdsettings.reader.ranging.mfcw_compsel.frs    =  NaN; % Hz reader sampling frequency
%     estimation window
stdsettings.reader.ranging.window = [0.25, 0.8]; % stationary parts of added MFCW signal 
%     distance estimation
stdsettings.reader.ranging.mfcw_calcdist.nc    = NaN; % number of secondary carriers
stdsettings.reader.ranging.mfcw_calcdist.fi    =  []; % Hz secondary carrier offset frequencies
stdsettings.reader.ranging.mfcw_calcdist.c_ord =  {}; % order in which to combine the carriers; empty: full permutation
stdsettings.reader.ranging.mfcw_calcdist.c     = stdsettings.c; % m/s speed of light


% *******************************************************************************************************
% channel reader <-> tag (only large-scale; no noise, no multipath propagation)

% large-scale and directivity
stdsettings.channel.largescale.type = 'log-dist'; % {'log-dist'}
stdsettings.channel.largescale.pl   =          2; % path-loss factor
stdsettings.channel.largescale.c    = stdsettings.c; % m/s speed of light
stdsettings.channel.largescale.f0   =         []; % Hz "center" frequencies for attenuation (assuming f0 >> fi)
stdsettings.channel.largescale.fs   =        NaN; % Hz sampling frequency
stdsettings.channel.directivity.txant  =      ''; % isotropic TX antenna
stdsettings.channel.directivity.txrot  =  [0, 0];
stdsettings.channel.directivity.rxant  =      ''; % isotropic RX antenna
stdsettings.channel.directivity.rxrot  =  [0, 0];
stdsettings.channel.directivity.dir_tx =  [0, 0];

% small-scale: OFF
stdsettings.channel.smallscale.on        = false; % switch smallscale model on/off
stdsettings.channel.smallscale.det       =   NaN; % purely deterministic?
stdsettings.channel.smallscale.seed      =   NaN; % random seed for non-deterministic setup (set to NaN for normal operation)
stdsettings.channel.smallscale.ensembles =   NaN; % number of independent realizations (set to 1!)
stdsettings.channel.smallscale.maxiter   =   NaN; % maximum number of iterations for PDP generation (sampling)
stdsettings.channel.smallscale.maxrays   =   NaN; % maximum number of paths (not including the LOS path)
stdsettings.channel.smallscale.k         =   NaN; % dB Ricean K-factor
stdsettings.channel.smallscale.trms      =   NaN; % s RMS delay spread
stdsettings.channel.smallscale.bw        =   NaN; % Hz bandwidth (w. oversampling)
stdsettings.channel.smallscale.fres      =   NaN; % initial frequency resolution in Hz (might be increased to meet .trms)
stdsettings.channel.smallscale.eps_k     =   NaN; % relative tolerable error in K-factor (linear)
stdsettings.channel.smallscale.eps_trms  =   NaN; % relative tolerable error in RMS delay spread
stdsettings.channel.smallscale.fs        =   NaN; % Hz sampling frequency

% noise: OFF
stdsettings.channel.noise.type =      'off'; % {'wgn', otherwise 'off'}
stdsettings.channel.noise.n0   = [NaN, NaN]; % single-sided noise density [dBm, per ? Hz]
stdsettings.channel.noise.fs   =        NaN; % Hz sampling frequency
stdsettings.channel.noise.frxs =        NaN; % Hz receiver (tag,reader) sampling frequency (for noise density -> power)

% Ricean K-Factor and RMS delay spread vs. distance
stdsettings.channel.v_dist = linspace(        0,       10, 100); % m (start at 0 m!)
stdsettings.channel.v_k    = logspace(log10(30), log10(6), 100); % dB
stdsettings.channel.v_trms = linspace(     1e-9,    20e-9, 100); % s


% *******************************************************************************************************
% reader

% oscillator
stdsettings.reader.oscillator.fcenter =         NaN; % Hz center frequency 
stdsettings.reader.oscillator.fstddev =        10e3; % Hz frequency standard deviation (Gaussian) ... phase noise
stdsettings.reader.oscillator.astddev =        1e-6; % amplitude standard deviation (center = 1; Gaussian)
stdsettings.reader.oscillator.snr     =         100; % SNR in dBc (to carrier)
stdsettings.reader.oscillator.length  =         NaN; % s length of carrier wave signal (will be defined later)
stdsettings.reader.oscillator.mode    =       'cos'; % cosine or sine wave
stdsettings.reader.oscillator.fs      =         NaN; % Hz sampling frequency

% modulation
stdsettings.reader.modulation.forcesettings =       true; % if true, no consistency/conformity checks are done SET TO TRUE!
stdsettings.reader.modulation.modulation    =   'PR-ASK'; % {'DSB-ASK', 'SSB-ASK', 'PR-ASK'}
stdsettings.reader.modulation.leadin        = 'preamble'; % {'preamble', 'frame-sync' or 'framesync'}
stdsettings.reader.modulation.tari          =      25e-6; % s duration of data-0 in (preferred: {6.25, 12.5, 25} us)
stdsettings.reader.modulation.lf            =        NaN; % Hz backscatter link frequency in kHz {40...640}
stdsettings.reader.modulation.dr            =          8; % divide ratio DR {64/3, 8}
stdsettings.reader.modulation.delimiter     =    12.5e-6; % (initial) delimiter duration in s ; -1: random
stdsettings.reader.modulation.x             =          1; % times tari {0.5...1 ; -1 means random} (data-1=tari+x)
stdsettings.reader.modulation.pw            =        0.5; % times tari {0.265 ... 0.525 ; -1 means random}
stdsettings.reader.modulation.moddepth      =         90; % percent {80...100 ; -1 means random}
stdsettings.reader.modulation.trf           =       0.25; % times tari {0...0.33 ; -1 means random}
stdsettings.reader.modulation.t0            =       1e-6; % s time delay (unmodulated carrier before and after modulation)
stdsettings.reader.modulation.fs            =        NaN; % Hz sampling frequency

% reader ouput stage (transmitter)
stdsettings.reader.transmitter.ptx     =             NaN; % W EIRP transmit power
stdsettings.reader.transmitter.en_nl   =           false; % enable nonlinear poweramp?
stdsettings.reader.transmitter.nl_char = 'readerchar_transmitter'; % filename for nonlinear power amplifier characteristic
stdsettings.reader.transmitter.nl_lim  =            0.95; % signal amplitude limits 0<=x<=1 (1: maximum limit of pwramp)
stdsettings.reader.transmitter.en_bp   =           false; %  enable bandpass true/false
stdsettings.reader.transmitter.bp_ord  =               4; % bandpass filter order
stdsettings.reader.transmitter.bp_att  =              80; % bandpass stopband attenuation in dB
stdsettings.reader.transmitter.bp_fcut = [700e6, 1200e6]; % Hz bandpass stopband edge frequencies in Hz
stdsettings.reader.transmitter.fs      =             NaN; % sampling frequency in Hz

% reader commands
%  query
stdsettings.reader.command.sel     =  'All'; % {'all', '~sl', 'sl'}
stdsettings.reader.command.session =   'S0'; % {'s0', 's1', 's2', 's3'}
stdsettings.reader.command.target  =    'A'; % {'a', 'b'}
stdsettings.reader.command.q       =      0; % {0, 1, 2, ..., 14, 15} 
stdsettings.reader.command.m       =      1; % {1, 2, 4, 8} (1:FM0, 2,4,8:Miller)};
stdsettings.reader.command.trext   =      0; % {0, 1} (pilot tone on/off)
stdsettings.reader.command.dr      = stdsettings.reader.modulation.dr; % divide ratio DR {64/3, 8}
%  ack
stdsettings.reader.command.rn16    = '0000'; % RN16 number (hex)

% reader input stage (receiver)
stdsettings.reader.receiver.en_bp  =           false; % enable bandpass filter?
stdsettings.reader.receiver.iirord =               4; % bandpasses filter order
stdsettings.reader.receiver.att    =              80; % stopband attenuation of bandpasses in dB
stdsettings.reader.receiver.fcut   = [700e6, 1200e6]; % Hz bandpass stopband edge frequencies
stdsettings.reader.receiver.fs     =             NaN; % Hz sampling frequency in MHz

% IQ demodulation
stdsettings.reader.demodulation.iirord     =           8; % anti-aliasing filter order
stdsettings.reader.demodulation.att        =         100; % dB stopband attenuation of lowpasses
stdsettings.reader.demodulation.frs        =         NaN; % Hz reader sampling frequency 
stdsettings.reader.demodulation.fcut       =         NaN; % Hz lowpass stopband edge frequency 
stdsettings.reader.demodulation.fs         =         NaN; % Hz sampling frequency
stdsettings.reader.demodulation.q          =         NaN; % bits quantization (set to NaN to deactivate)
stdsettings.reader.demodulation.len_rs     =         NaN; % length of output signal 

% analog path frequency response estimation
stdsettings.reader.analogpathest.f0   =         NaN; % Hz main carrier frequency 
stdsettings.reader.analogpathest.fs   =         NaN; % Hz sampling frequency
%     error in estimates: ~Gaussian [mult. factor (syst. err.), std. dev. factor (rand err.)]
stdsettings.reader.analogpathest.err  =      [1, 0]; % no error: [1, 0]


% *******************************************************************************************************
% tag

% tag identification
stdsettings.tag.id.rn16 = '208C'; % RN16 number (hex)
stdsettings.tag.id.pc   = '0006'; % protocol-control (PC) bits (hex)
stdsettings.tag.id.epc  = '0013 1F4A 0000 0000 0000 0017'; % electronic product code (EPC) (hex)

% tag state (may / may not be identical to EPCglobal)
stdsettings.tag.state.powered =  0; % tag is powered?
stdsettings.tag.state.epc     = ''; % state of tag according to EPCglobal Class-1 Gen-2

% clock
stdsettings.tag.clock.fcenter   =         3e6; % Hz center ("clock") frequency in MHz ... set relatively high
stdsettings.tag.clock.fsigma    =           0; % percent sampling time variation (sigma)
stdsettings.tag.clock.phi0      =          -1; % deg phase shift (0<=phi0<360) ; -1: random
stdsettings.tag.clock.length    =         NaN; % samples maximum sample time 
stdsettings.tag.clock.mode      =      'udef'; % {'fs', 'fclk'} 
stdsettings.tag.clock.fs        =         NaN; % Hz sampling frequency

% modulation
stdsettings.tag.modulation.force       =  true; % force modulation / filtering (if tag is not functional)
stdsettings.tag.modulation.ptag_f      =   NaN; % forced power level (NaN for normal operation)
stdsettings.tag.modulation.pfnc_bias   =  0.98; % bias towards "unmodulated" Pin and Pic 0...1 (0:mod+unmod, 1:unmod); MFCW: short-term
stdsettings.tag.modulation.t0          =  1e-4; % s initial delay
stdsettings.tag.modulation.lf          =   stdsettings.reader.modulation.lf; % Hz backscatter link frequency (will be set by reader cmd for return link)
stdsettings.tag.modulation.charfile    =    ''; % reflection coeff. mat-file 
stdsettings.tag.modulation.pwrchar     =      'Vdda_AVG'; % power supply characteristic
stdsettings.tag.modulation.pwrcharfile = 'tagchar_power'; % tag power-supply characteristics filename
stdsettings.tag.modulation.nfilt       =    32; % filterbank channel filter length
stdsettings.tag.modulation.nfft        =  4096; % filterbank channels
stdsettings.tag.modulation.fclk        = stdsettings.tag.clock.fcenter; % Hz tag clock frequency
stdsettings.tag.modulation.fc          = stdsettings.reader.oscillator.fcenter; % Hz carrier frequency
stdsettings.tag.modulation.fs          =   NaN;  % Hz sampling frequency

% demodulation
stdsettings.tag.demodulation.agc_sat        =           1; % keep functional for Pcarrier, Vdda out of range (AGC)
stdsettings.tag.demodulation.vdda           =         NaN; % V power supply (Vdda) voltage (scalar); Vgnd = Agnd = 0 V
stdsettings.tag.demodulation.tau            =   169.12e-9; % s time constant of input RC filter
stdsettings.tag.demodulation.force_exact_rc =       false; % force recalculation of RC filter coeff. if loaded ones are not an exact match
stdsettings.tag.demodulation.agc_len        =         250; % periods of fc window length for AGC power detection (fcn. EST_POWER)
stdsettings.tag.demodulation.agc_ol         =          50; % percent overlapping between AGC power detection windows (fcn. EST_POWER)
stdsettings.tag.demodulation.agc_nbins      =         100; % number of bins for AGC amplitude histogram
stdsettings.tag.demodulation.dvpeak         =       58e-3; % V steady state venv - vpeak
stdsettings.tag.demodulation.slewrate       =      .045e6; % V/s maximum slew rate for vpeak
stdsettings.tag.demodulation.ar             =        5e-5; % V rise factor for vpeak (empirical)
stdsettings.tag.demodulation.af             = stdsettings.tag.demodulation.ar*1e-4; % V fall factor for vpeak (empirical)
stdsettings.tag.demodulation.debouncelength =           2; % length of debounce filter (MA filter)
stdsettings.tag.demodulation.fs             =         NaN; % Hz sampling frequency
stdsettings.tag.demodulation.fc             = stdsettings.reader.oscillator.fcenter; % Hz carrier (center) frequency
stdsettings.tag.demodulation.fclk           = stdsettings.tag.clock.fcenter; % Hz tag clock (center) frequency
%     downsampling for performance reasons (use with care!)
stdsettings.tag.demodulation.turbo          =     true; % downsampling prior to vpeak calculation if true
stdsettings.tag.demodulation.turbo_osfactor =      100; % oversampling rate (critical sampling: 2*fclk) for "turbo" downsampling

% encoding
stdsettings.tag.encoding.trext = NaN; % pilot tone on/off
stdsettings.tag.encoding.m     = NaN; % (1:FM0, 2,4,8:Miller)

% decoding
stdsettings.tag.decoding.fclk = stdsettings.tag.clock.fcenter; % Hz center ("clock") frequency in MHz


% *******************************************************************************************************
% *******************************************************************************************************
% GENERATE SETTINGS
disp('Generating settings.');


% *******************************************************************************************************
% complete/check partitioning 

% partitioning
if length(partitioning.names) ~= length(partitioning.runs)
   error('Length of entries of struct partitioning do not match. Check partitioning.');
end
%     indices [0, end-of-block]
partitioning.indices = [0, cumsum(partitioning.runs)]; 


% *******************************************************************************************************
% basic (random) settings
   
for p = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(p) : partitioning.indices(p+1)
      % general
      settings{j}.rand.fs = randsettings.fs(1) + rand*(randsettings.fs(2) - randsettings.fs(1));
      settings{j}.c       = stdsettings.c;
      % reader
      settings{j}.rand.fc    = randsettings.fc (1) + rand*(randsettings.fc (2) - randsettings.fc (1));
      settings{j}.rand.frs   = randsettings.frs(1) + rand*(randsettings.frs(2) - randsettings.frs(1));
      settings{j}.rand.ptx   = randsettings.ptx(1) + rand*(randsettings.ptx(2) - randsettings.ptx(1));
      settings{j}.rand.pos_r = 0; % set reader position to zero for simplicity
      settings{j}.rand.qy_s   = randsettings.qy_s { intrand([1, length(randsettings.qy_s )]) };
      settings{j}.rand.qy_t   = randsettings.qy_t { intrand([1, length(randsettings.qy_t )]) };
      settings{j}.rand.qy_q   = intrand(randsettings.qy_q);
      %     make sure oversampling rate ist not close to an integer 
      %     (would result in very low beat frequency which cannot be handled by tag_demodulation)
      settings{j}.rand.fs = settings{j}.rand.fs + ...
         settings{j}.rand.fc/2 - mod(settings{j}.rand.fs, settings{j}.rand.fc);
      % ranging (part of reader)
      settings{j}.rand.mfcw.nc   = intrand(randsettings.nc);
      settings{j}.rand.mfcw.fi   = cumsum( randsettings.fi(1) +...
         rand(1,settings{j}.rand.mfcw.nc)*(randsettings.fi(2) - randsettings.fi(1)) );
      settings{j}.rand.mfcw.vari = randsettings.vari(1) +...
         rand(1,settings{j}.rand.mfcw.nc)*(randsettings.vari(2) - randsettings.vari(1));
      % tag
      settings{j}.rand.fm      = randsettings.fm  (1) + rand*(randsettings.fm  (2) - randsettings.fm  (1));
      settings{j}.rand.pwrchar = randsettings.pwrc{ intrand([1, length(randsettings.pwrc)]) };
      settings{j}.rand.modchar = randsettings.modc{ intrand([1, length(randsettings.modc)]) };
      if strcmpi(partitioning.names{p}, 'MFCW')
         settings{j}.rand.rn16  = ''; % zero-length RN16 to speed up the tests / limit RAM usage
         settings{j}.rand.pos_t = randsettings.dist(1) + rand*(randsettings.dist(2) - randsettings.dist(1));
      else
         settings{j}.rand.rn16  = dec2hex(round(rand*2^16), 4);
         settings{j}.rand.pos_t = 1; % fixed 1m distance for all protocol tests
      end
      % round distance to sampling interval
      settings{j}.rand.pos_t = round(settings{j}.rand.pos_t/settings{j}.c*settings{j}.rand.fs)*settings{j}.c/settings{j}.rand.fs;
   end
end


% *******************************************************************************************************
% complete setup and create individuals (based upon sim_mfcw_ranging alpha 1.2)

for p = 1 : length(partitioning.names)
   for j = 1 + partitioning.indices(p) : partitioning.indices(p+1)
      % copy from standardsettings
      settings{j}.reader  = stdsettings.reader;
      
      % **************************************************
      % reader
      
      % component setup
      settings{j}.reader.oscillator.fcenter = settings{j}.rand.fc;
      settings{j}.reader.oscillator.fs      = settings{j}.rand.fs;
      settings{j}.reader.modulation.fs      = settings{j}.rand.fs;
      settings{j}.reader.modulation.lf      = settings{j}.rand.fm;
      settings{j}.reader.analogpathest.f0   = settings{j}.rand.fc;
      settings{j}.reader.analogpathest.fs   = settings{j}.rand.fs;
      settings{j}.reader.transmitter.ptx    = settings{j}.rand.ptx;
      settings{j}.reader.transmitter.fs     = settings{j}.rand.fs;
      settings{j}.reader.receiver.fs        = settings{j}.rand.fs;
      settings{j}.reader.analogpathest.fs   = settings{j}.rand.fs;
      settings{j}.reader.demodulation.frs   = settings{j}.rand.frs;
      settings{j}.reader.demodulation.fs    = settings{j}.rand.fs;
      settings{j}.reader.demodulation.fcut  =...
         3 * ( max(abs(settings{j}.rand.mfcw.fi)) + settings{j}.reader.modulation.lf );
      
      % MFCW ranging setup (oscillator setup will be done later)
      %     general
      fm = settings{j}.reader.modulation.lf;
      fi = settings{j}.rand.mfcw.fi;
      settings{j}.rand.mfcw.fi_rs = subsampling(fi, settings{j}.reader.demodulation.frs);
      if any(abs(diff(abs(settings{j}.rand.mfcw.fi_rs))) <= 3*fm)
         error('Detected overlapping of carriers/sidebands. Check frequency setup (subsampling?).')
      end
      settings{j}.reader.ranging.freq = [-fm, 0, fm, fi-fm, fi, fi+fm]; % for analog input stage estimator
      %     secondary carrier generation
      settings{j}.reader.ranging.mfcw_addseccarriers.f0   = settings{j}.rand.fc;
      settings{j}.reader.ranging.mfcw_addseccarriers.nc   = settings{j}.rand.mfcw.nc;
      settings{j}.reader.ranging.mfcw_addseccarriers.fi   = settings{j}.rand.mfcw.fi;
      settings{j}.reader.ranging.mfcw_addseccarriers.vari = settings{j}.rand.mfcw.vari;
      settings{j}.reader.ranging.mfcw_addseccarriers.fs   = settings{j}.rand.fs;
      %     component selection
      settings{j}.reader.ranging.mfcw_compsel.nc  = settings{j}.rand.mfcw.nc;
      settings{j}.reader.ranging.mfcw_compsel.fi  = settings{j}.rand.mfcw.fi_rs;
      settings{j}.reader.ranging.mfcw_compsel.lf  = fm;
      settings{j}.reader.ranging.mfcw_compsel.bw  = fm/2;
      settings{j}.reader.ranging.mfcw_compsel.frs = settings{j}.reader.demodulation.frs;
      %     distance estimator setup
      settings{j}.reader.ranging.mfcw_calcdist.nc = settings{j}.rand.mfcw.nc;
      settings{j}.reader.ranging.mfcw_calcdist.fi = settings{j}.rand.mfcw.fi;
      %     cleanup
      clear('fm', 'fi');
      
      % test specific settings
      if strcmpi(partitioning.names{p}, 'MFCW') % long pilot tone
         settings{j}.reader.command.trext = 1; % pilot tone on
         settings{j}.reader.command.m     = 8; % long pilot tone
      else
         settings{j}.reader.command.trext = 0; % no pilot tone and
         settings{j}.reader.command.m     = 1; % FM0 encoding (speed)
      end
      
      % mfcw_theorysim (for MFCW tests)
      if strcmpi(partitioning.names{p}, 'MFCW')
         % general
         settings{j}.reader.mfcw_theorysim.silent  =  true;
         settings{j}.reader.mfcw_theorysim.verbose = false;
         % frequency and time setup
         settings{j}.reader.mfcw_theorysim.f_s   = settings{j}.rand.fs;
         settings{j}.reader.mfcw_theorysim.f_rs  = settings{j}.reader.demodulation.frs;
         settings{j}.reader.mfcw_theorysim.t_max = NaN; % set by test
         % carrier setup<
         settings{j}.reader.mfcw_theorysim.carriers.N_c  = settings{j}.rand.mfcw.nc + 1;
         A_i  = [1, settings{j}.rand.mfcw.vari];
         settings{j}.reader.mfcw_theorysim.carriers.A_i  = sqrt(2*settings{j}.rand.ptx) * A_i / sqrt(sum(A_i.^2)); % normalized output power Ptx
         settings{j}.reader.mfcw_theorysim.carriers.f_i  = [0, settings{j}.rand.mfcw.fi];
         settings{j}.reader.mfcw_theorysim.carriers.f_c  = settings{j}.rand.fc;
         settings{j}.reader.mfcw_theorysim.carriers.comb = 'sim';
         % tag setup
         settings{j}.reader.mfcw_theorysim.tag.f_m         = settings{j}.reader.modulation.lf;
         settings{j}.reader.mfcw_theorysim.tag.tagchar_mod = '';% set by test;
         settings{j}.reader.mfcw_theorysim.tag.p_av        = NaN;
         % channel setup
         settings{j}.reader.mfcw_theorysim.dist          = NaN; % set by test
         settings{j}.reader.mfcw_theorysim.awgn.N_0      = [-Inf, Inf]; % no noise
         settings{j}.reader.mfcw_theorysim.awgn.corr_ci  = 0;
         settings{j}.reader.mfcw_theorysim.awgn.corr_iq  = 0;
         settings{j}.reader.mfcw_theorysim.feedback.gain = 0; % no feedback
         % filter setups
         settings{j}.reader.mfcw_theorysim.filter.ord_iq  = settings{j}.reader.demodulation.iirord;
         settings{j}.reader.mfcw_theorysim.filter.att_iq  = settings{j}.reader.demodulation.att;
         settings{j}.reader.mfcw_theorysim.filter.ord_sel = settings{j}.reader.ranging.mfcw_compsel.iirord;
         settings{j}.reader.mfcw_theorysim.filter.att_sel = settings{j}.reader.ranging.mfcw_compsel.att;
         % system estimation
         settings{j}.reader.mfcw_theorysim.sysest.A_err = [1, 0]; % no errors
         settings{j}.reader.mfcw_theorysim.sysest.G_err = [1, 0];
         % random parameters [mean, std]
         settings{j}.reader.mfcw_theorysim.random.arg_A_i = [0,0];
         settings{j}.reader.mfcw_theorysim.random.ms_del  = [0,0];
         settings{j}.reader.mfcw_theorysim.random.phi_m   = [0,0];
         settings{j}.reader.mfcw_theorysim.random.fb_del  = [0,0];
         % other
         settings{j}.reader.mfcw_theorysim.c = settings{j}.c;
      end
      
      
      % **************************************************
      % tag
      
      % copy from standardsettings
      settings{j}.tag = stdsettings.tag;
      
      % tag component setup
      settings{j}.tag.id.rn16             = settings{j}.rand.rn16;
      settings{j}.tag.clock.fs            = settings{j}.rand.fs;
      settings{j}.tag.modulation.fs       = settings{j}.rand.fs;
      settings{j}.tag.demodulation.fs     = settings{j}.rand.fs;
      settings{j}.tag.modulation.fc       = settings{j}.rand.fc;
      settings{j}.tag.modulation.charfile = settings{j}.rand.modchar;
      settings{j}.tag.modulation.pwrchar  = settings{j}.rand.pwrchar;
      settings{j}.tag.demodulation.fc     = settings{j}.rand.fc;
      settings{j}.tag.modulation.lf       = settings{j}.reader.modulation.lf;
      
      % make tag clock frequency a multiple of LF
      settings{j}.tag.clock.fcenter = settings{j}.tag.clock.fcenter +...
         settings{j}.rand.fm - mod(settings{j}.tag.clock.fcenter, settings{j}.rand.fm);
      
      % **************************************************
      % channels reader <-> tag
      
      % copy from standardsettings
      %     reader -> tag
      settings{j}.channel_rt{1,1}.largescale  = stdsettings.channel.largescale;
      settings{j}.channel_rt{1,1}.directivity = stdsettings.channel.directivity;
      settings{j}.channel_rt{1,1}.smallscale  = stdsettings.channel.smallscale;
      settings{j}.channel_rt{1,1}.noise       = stdsettings.channel.noise;
      %     tag -> reader
      settings{j}.channel_tr{1,1}.largescale  = stdsettings.channel.largescale;
      settings{j}.channel_tr{1,1}.directivity = stdsettings.channel.directivity;
      settings{j}.channel_tr{1,1}.smallscale  = stdsettings.channel.smallscale;
      settings{j}.channel_tr{1,1}.noise       = stdsettings.channel.noise;
      % largescale
      settings{j}.channel_rt{1,1}.largescale.f0 = settings{j}.rand.fc;
      settings{j}.channel_tr{1,1}.largescale.f0 = settings{j}.rand.fc;
      % largescale
      settings{j}.channel_rt{1,1}.largescale.fs = settings{j}.rand.fs;
      settings{j}.channel_tr{1,1}.largescale.fs = settings{j}.rand.fs;
      % smallscale
      settings{j}.channel_rt{1,1}.smallscale.bw   = 3 * ( 2*settings{j}.reader.modulation.lf +...
         max([0, settings{j}.rand.mfcw.fi]) - min([0, settings{j}.rand.mfcw.fi]) );
      settings{j}.channel_tr{1,1}.smallscale.bw   = 3 * ( 2*settings{j}.reader.modulation.lf +...
         max([0, settings{j}.rand.mfcw.fi]) - min([0, settings{j}.rand.mfcw.fi]) );
      settings{j}.channel_rt{1,1}.smallscale.fres = min(diff(sort([0, settings{j}.rand.mfcw.fi]))); % smallest carrier spacing
      settings{j}.channel_tr{1,1}.smallscale.fres = min(diff(sort([0, settings{j}.rand.mfcw.fi])));
      % noise
      settings{j}.channel_rt{1,1}.noise.frxs = settings{j}.tag.clock.fcenter; % receiver is tag
      settings{j}.channel_tr{1,1}.noise.frxs = settings{j}.reader.demodulation.frs; % receiver is reader
      
      % add all distance-dependent settings
      settings{j}.channel_rt = channel_newpos(settings{j}.channel_rt, ...
         settings{j}.rand.pos_r, settings{j}.rand.pos_t, stdsettings.channel);
      settings{j}.channel_tr = channel_newpos(settings{j}.channel_tr, ...
         settings{j}.rand.pos_t, settings{j}.rand.pos_r, stdsettings.channel);
      
      % check for possible phase ambiguities (on the safe side)
      if  5 * settings{j}.channel_rt{1,1}.largescale.dist >= settings{j}.c / max(diff([0, settings{j}.rand.mfcw.fi]))
         error('Possible phase ambiguities. Check setup (c, dist, fi).')
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
matfilename    = filename; 
characteristic = 'settings and results for selftest: reader_main and tag_main'; 
createdon      = datestr(now, 0); 
createdby      = sprintf('%s.m, rev.: %s', mfilename, version);
usedmatfiles   = ''; 

% save (overwrite if exists)
save(fullfile(directory, filename), 'matfilename','characteristic','createdby','usedmatfiles',...
               'settings', 'partitioning');
disp(sprintf('\nSettings and results have been saved in %s/%s.mat', pwd, fullfile(directory, filename) ));

