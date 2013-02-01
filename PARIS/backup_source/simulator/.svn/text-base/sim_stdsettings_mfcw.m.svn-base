% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% standard setup for mfcw ranging (for later customization)
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
% version = sim_stdsettings()
%    Just returns the version number (string).
% stdsettings = sim_stdsettings(settings)
%    Returns a standard setup for customization. Note that the returned setup is incomplete; undefined
%    values are represented by NaNs or empty arrays.
%
%
% ***** Interface definition *****
% function stdsettings = sim_stdsettings(settings)
%    settings    struct containing basic settings
%       .fs         sampling frequency in Hz
%       .c          speed of light in m/s
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

function stdsettings = sim_stdsettings_mfcw(settings)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   stdsettings = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% *******************************************************************************************************
% SETTINGS

% *******************************************************************************************************
% ranging

% assumed tolerance for preamble position
% stdsettings.ranging.preamble_tol = max(cell2mat(tagpool.t0)); % us

% multifrequency-CW-radar (mfcw)
%     secondary carrier generation
stdsettings.reader.ranging.mfcw_addseccarriers.nc   =         NaN; % number of secondary carriers
stdsettings.reader.ranging.mfcw_addseccarriers.fi   =          []; % Hz secondary carrier frequencies
stdsettings.reader.ranging.mfcw_addseccarriers.vari =          []; % secondary carrier variances
stdsettings.reader.ranging.mfcw_addseccarriers.f0   =         NaN; % Hz primary carrier frequency
stdsettings.reader.ranging.mfcw_addseccarriers.fs   = settings.fs; % Hz sampling frequency
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
stdsettings.reader.ranging.mfcw_calcdist.c_ord =  {}; % order in which to combine the carriers {2} for 2FCW, ...
stdsettings.reader.ranging.mfcw_calcdist.c     = settings.c; % m/s speed of light


% *******************************************************************************************************
% channel reader <-> tag

% distance reader -> tag
%              tag1, tag2, ...
%     reader1:
%     reader2:
%     reader3:
%        :
stdsettings.dist_rt = []; % m
%     round to sampling resolution
stdsettings.dist_rt = round(stdsettings.dist_rt/settings.c*settings.fs)*settings.c/settings.fs;

% distance tag -> reader
%           reader1, reader2, ...
%     tag1:
%     tag2:
%     tag3:
stdsettings.dist_tr = stdsettings.dist_rt';

% Ricean K-Factor and RMS delay spread vs. distance
stdsettings.channel.v_dist = linspace(        0,       10, 1000); % m (start at 0 m!)
stdsettings.channel.v_k    = logspace(log10(30), log10(6), 1000); % dB
stdsettings.channel.v_trms = linspace(     1e-9,    20e-9, 1000); % s

% channel type and basic settings
stdsettings.channel.probemode =    false; % only record channels if true; do not calculate output
stdsettings.channel.type =        'room'; % channel type (does directivity apply to scattered?) {'room' -> no, 'outdoor' -> yes}) 
stdsettings.channel.c    =    settings.c; % m/s speed of light
stdsettings.channel.fs   =   settings.fs; % Hz sampling frequency

% large-scale
%     log-dist:  log-distance model
%     logn-shad: log-distance with log-normal shadowing ... NOT IMPLEMENTED YET
stdsettings.channel.largescale.type =    'log-dist'; % {'log-dist'}
stdsettings.channel.largescale.pl   =             2; % path-loss factor
stdsettings.channel.largescale.c    =    settings.c; % m/s speed of light
stdsettings.channel.largescale.f0   =            []; % Hz "center" frequencies for attenuation (assuming f0 >> fi)
stdsettings.channel.largescale.fs   =   settings.fs; % Hz sampling frequency

% directivity
stdsettings.channel.directivity.txant  =     ''; % filename of TX antenna characteristic (isotropic if empty)
stdsettings.channel.directivity.txrot  = [0, 0]; % rotation of TX antenna [azimuth, elevation] in degree ([0,0]: maximum in positive x-axis)
stdsettings.channel.directivity.rxant  =     ''; % filename of RX antenna characteristic (isotropic if empty)
stdsettings.channel.directivity.rxrot  = [0, 0]; % rotation of RX antenna [azimuth, elevation] in degree ([0,0]: maximum in positive x-axis)
stdsettings.channel.directivity.dir_tx =     []; % direction of channel [azimuth, elevation] in degree (from TX antenna)

% small-scale
stdsettings.channel.smallscale.on        =     1; % switch smallscale model on/off
stdsettings.channel.smallscale.det       =     0; % purely deterministic?
stdsettings.channel.smallscale.seed      =   NaN; % random seed for non-deterministic setup (set to NaN for normal operation)
stdsettings.channel.smallscale.ensembles =     1; % number of independent realizations (set to 1!)
stdsettings.channel.smallscale.maxiter   =  1000; % maximum number of iterations for PDP generation (sampling)
stdsettings.channel.smallscale.maxrays   =  1000; % maximum number of paths (not including the LOS path)
stdsettings.channel.smallscale.k         =   NaN; % dB Ricean K-factor
stdsettings.channel.smallscale.trms      =   NaN; % s RMS delay spread
stdsettings.channel.smallscale.bw        =   NaN; % Hz bandwidth (w. oversampling)
stdsettings.channel.smallscale.fres      =   NaN; % initial frequency resolution in Hz (might be increased to meet .trms)
stdsettings.channel.smallscale.eps_k     =  1e-3; % tolerable relative error of K-factor (linear)
stdsettings.channel.smallscale.eps_trms  =  1e-2; % tolerable relative error of RMS delay spread
stdsettings.channel.smallscale.fs        = settings.fs; % Hz sampling frequency

% noise
%     off: no noise
%     wgn: white gaussian noise
stdsettings.channel.noise.type =       'off'; % {'wgn', otherwise 'off'}
stdsettings.channel.noise.n0   = [-47, 3000]; % single-sided noise density [dBm, per ? Hz]
stdsettings.channel.noise.fs   = settings.fs; % Hz sampling frequency
stdsettings.channel.noise.frxs =         NaN; % Hz receiver (tag,reader) sampling frequency (for noise density -> power)


% *******************************************************************************************************
% channel reader <-> reader (feedback)

% start from channel setup
stdsettings.feedback = stdsettings.channel;

% distance reader <-> reader
%                 rx reader1, rx reader2, ...
%     tx reader1:
%     tx reader2:
%     tx reader4:
%        :
stdsettings.dist_rr = []; % m (main diagonal should contain zeros ... handled below)
%     round to sampling resolution
stdsettings.dist_rr = round(stdsettings.dist_rr/settings.c*settings.fs)*settings.c/settings.fs;

% gain/delay for direct feedback
stdsettings.feedback.fb_att =   50; % dB direct feedback attenuation (active carrier suppression)
stdsettings.feedback.fb_del = 3e-9; % s direct feedback delay (will be rounded to sampling time)

% feedback adds noise (see signal model)
stdsettings.feedback.noise.type = 'wgn';


% *******************************************************************************************************
% reader

% oscillator
stdsettings.reader.oscillator.fcenter =         NaN; % Hz center frequency 
stdsettings.reader.oscillator.fstddev =           0; % Hz frequency standard deviation (Gaussian) ... phase noise
stdsettings.reader.oscillator.astddev =           0; % amplitude standard deviation (center = 1; Gaussian)
stdsettings.reader.oscillator.snr     =         Inf; % SNR in dBc (to carrier)
stdsettings.reader.oscillator.length  =         NaN; % s length of carrier wave signal (will be defined later)
stdsettings.reader.oscillator.mode    =       'cos'; % cosine or sine wave
stdsettings.reader.oscillator.fs      = settings.fs; % Hz sampling frequency

% modulation
stdsettings.reader.modulation.forcesettings =        false; % if true, no consistency/conformity checks are done
stdsettings.reader.modulation.modulation    =     'PR-ASK'; % {'DSB-ASK', 'SSB-ASK', 'PR-ASK'}
stdsettings.reader.modulation.leadin        =   'preamble'; % {'preamble', 'frame-sync' or 'framesync'}
stdsettings.reader.modulation.tari          =        25e-6; % s duration of data-0 in (preferred: {6.25, 12.5, 25} us)
stdsettings.reader.modulation.lf            =         48e3; % Hz backscatter link frequency in kHz {40...640}
stdsettings.reader.modulation.dr            =            8; % divide ratio DR {64/3, 8}
stdsettings.reader.modulation.delimiter     =      12.5e-6; % (initial) delimiter duration in s ; -1: random
stdsettings.reader.modulation.x             =            1; % times tari {0.5...1 ; -1 means random} (data-1=tari+x)
stdsettings.reader.modulation.pw            =          0.5; % times tari {0.265 ... 0.525 ; -1 means random}
stdsettings.reader.modulation.moddepth      =           90; % percent {80...100 ; -1 means random}
stdsettings.reader.modulation.trf           =         0.25; % times tari {0...0.33 ; -1 means random}
stdsettings.reader.modulation.t0            =         1e-6; % s time delay (unmodulated carrier before and after modulation)
stdsettings.reader.modulation.fs            =  settings.fs; % Hz sampling frequency

% reader ouput stage (transmitter)
stdsettings.reader.transmitter.ptx     =             NaN; % W EIRP transmit power
stdsettings.reader.transmitter.en_nl   =           false; % enable nonlinear poweramp?
stdsettings.reader.transmitter.nl_char = 'readerchar_transmitter'; % filename for nonlinear power amplifier characteristic
stdsettings.reader.transmitter.nl_lim  =            0.95; % signal amplitude limits 0<=x<=1 (1: maximum limit of pwramp)
stdsettings.reader.transmitter.en_bp   =           false; %  enable bandpass true/false
stdsettings.reader.transmitter.bp_ord  =               4; % bandpass filter order
stdsettings.reader.transmitter.bp_att  =              80; % bandpass stopband attenuation in dB
stdsettings.reader.transmitter.bp_fcut = [700e6, 1200e6]; % Hz bandpass stopband edge frequencies in Hz
stdsettings.reader.transmitter.fs      =     settings.fs; % sampling frequency in Hz

% reader commands
%  query
stdsettings.reader.command.sel     =  'All'; % {'all', '~sl', 'sl'}
stdsettings.reader.command.session =   'S0'; % {'s0', 's1', 's2', 's3'}
stdsettings.reader.command.target  =    'A'; % {'a', 'b'}
stdsettings.reader.command.q       =      0; % {0, 1, 2, ..., 14, 15} 
stdsettings.reader.command.m       =      8; % {1, 2, 4, 8} (1:FM0, 2,4,8:Miller)};
stdsettings.reader.command.trext   =      1; % {0, 1} (pilot tone on/off)
stdsettings.reader.command.dr      = stdsettings.reader.modulation.dr; % divide ratio DR {64/3, 8}
%  ack
stdsettings.reader.command.rn16    = '0000'; % RN16 number (hex)

% reader input stage (receiver)
stdsettings.reader.receiver.en_bp  =           false; % enable bandpass filter?
stdsettings.reader.receiver.iirord =               4; % bandpasses filter order
stdsettings.reader.receiver.att    =              80; % stopband attenuation of bandpasses in dB
stdsettings.reader.receiver.fcut   = [700e6, 1200e6]; % Hz bandpass stopband edge frequencies
stdsettings.reader.receiver.fs     =     settings.fs; % Hz sampling frequency in MHz

% IQ demodulation
stdsettings.reader.demodulation.iirord     =           8; % anti-aliasing filter order
stdsettings.reader.demodulation.att        =         100; % dB stopband attenuation of lowpasses
stdsettings.reader.demodulation.frs        =         NaN; % Hz reader sampling frequency 
stdsettings.reader.demodulation.fcut       =         NaN; % Hz lowpass stopband edge frequency 
stdsettings.reader.demodulation.fs         = settings.fs; % Hz sampling frequency
stdsettings.reader.demodulation.q          =         NaN; % bits quantization (set to NaN to deactivate)
stdsettings.reader.demodulation.len_rs     =         NaN; % length of output signal 

% analog path frequency response estimation
stdsettings.reader.analogpathest.f0   =         NaN; % Hz main carrier frequency 
stdsettings.reader.analogpathest.fs   = settings.fs; % Hz sampling frequency
%     error in estimates: ~Gaussian [mult. factor (syst. err.), std. dev. factor (rand err.)]
stdsettings.reader.analogpathest.err  = [1, 0]; % no error: [1, 0]


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
stdsettings.tag.clock.fcenter   =      1.92e6; % Hz center ("clock") frequency in MHz
stdsettings.tag.clock.fsigma    =         5/3; % percent sampling time variation (sigma)
stdsettings.tag.clock.phi0      =          -1; % deg phase shift (0<=phi0<360) ; -1: random
stdsettings.tag.clock.length    =         NaN; % samples maximum sample time 
stdsettings.tag.clock.mode      =      'udef'; % {'fs', 'fclk'} 
stdsettings.tag.clock.fs        = settings.fs; % Hz sampling frequency

% modulation
stdsettings.tag.modulation.force       =  true; % force modulation / filtering (if tag is not functional)
stdsettings.tag.modulation.ptag_f      =   NaN; % forced power level (NaN for normal operation)
stdsettings.tag.modulation.pfnc_bias   =  0.98; % bias towards "unmodulated" Pin and Pic 0...1 (0:mod+unmod, 1:unmod); MFCW: short-term
stdsettings.tag.modulation.t0          =     0; % s initial delay
stdsettings.tag.modulation.lf          =   stdsettings.reader.modulation.lf; % Hz backscatter link frequency (will be set by reader cmd for return link)
stdsettings.tag.modulation.charfile    =    ''; % reflection coeff. mat-file 
stdsettings.tag.modulation.pwrchar     =      'Vdda_AVG'; % power supply characteristic
stdsettings.tag.modulation.pwrcharfile = 'tagchar_power'; % tag power-supply characteristics filename
stdsettings.tag.modulation.nfilt       =    32; % filterbank channel filter length
stdsettings.tag.modulation.nfft        =  4096; % filterbank channels
stdsettings.tag.modulation.fclk        = stdsettings.tag.clock.fcenter; % Hz tag clock frequency
stdsettings.tag.modulation.fc          = stdsettings.reader.oscillator.fcenter; % Hz carrier frequency
stdsettings.tag.modulation.fs          = settings.fs;  % Hz sampling frequency

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
stdsettings.tag.demodulation.fs             = settings.fs; % Hz sampling frequency
stdsettings.tag.demodulation.fc             = stdsettings.reader.oscillator.fcenter; % Hz carrier (center) frequency
stdsettings.tag.demodulation.fclk           = stdsettings.tag.clock.fcenter; % Hz tag clock (center) frequency
%     downsampling for performance reasons (use with care!)
stdsettings.tag.demodulation.turbo          =     true; % downsampling prior to vpeak calculation if true
stdsettings.tag.demodulation.turbo_osfactor =       25; % oversampling rate (critical sampling: 2) for "turbo" downsampling

% encoding
stdsettings.tag.encoding.trext = stdsettings.reader.command.trext; % pilot tone on/off
stdsettings.tag.encoding.m     = stdsettings.reader.command.m; % (1:FM0, 2,4,8:Miller)

% decoding
stdsettings.tag.decoding.fclk = stdsettings.tag.clock.fcenter; % Hz center ("clock") frequency in MHz
