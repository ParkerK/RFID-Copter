% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% reader - transmitter (output) stage (power amplifier, antenna + antenna gain)
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
% version = reader_transmitter()
%    Just returns the version number (string).
% Hd = reader_transmitter(settings)
%    Returns the overall (linear) filter object of this function (dfilt). 
%    ATTENTION: Note that this model does not include the nonlinear amplifier (purely linear model).
%               Furthermore this estimate assumes a normalized input signal power VAR(SIGNAL)=1!
% signal = reader_transmitter(signal, settings)
%    Amplifies and filters SIGNAL according to SETTINGS, emulating a linear/nonlinear (optional) power
%    amplifier with optional subsequent bandpass filtering. The function creates an exact output variance
%    for both (linear and nonlinear) amplifier types.
%
%
% ***** Interface definition *****
% function signal = reader_receiver(signal, settings)
%    signal        input signal
%    settings      struct containing settings
%       .ptx          desired transmit power (after bandpass for simplicity of setup)
%       .en_nl        pwramp: enable nonlinear power amplifier true/false (false: linear amplifier)
%       .nl_char      pwramp: filename for nonlinear power amplifier characteristic
%       .nl_lim       pwramp: symmetrical input signal amplitude limits 0<=x<=1 (1: max limits of pwramp)
%       .en_bp        bandpass: enable true/false
%       .bp_ord       bandpass: IIR filter order
%       .bp_att       bandpass: stopband attenuation in dB
%       .bp_fcut      bandpass: stopband edge frequencies in Hz
%       .fs           sampling frequency in Hz
%
%    signal   output signal according to SETTINGS
%             -> alternatively: filter object (overall dfilt filter structure of this function) 
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

function signal = reader_transmitter(varargin)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   signal = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% tolerance for symmetry of poweramplifier amplitude characteristic (has to be symmetrical around zero)
internalsettings.pwramp_symtol = 10 * eps;


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% number of input parameters
if nargin == 1
   justgauge = true; % just return linear model of this fcn
   settings  = varargin{1};
elseif nargin == 2
   justgauge   = false;
   signal      = varargin{1}(:);
   settings    = varargin{2};
else
   criterr('Wrong number of input arguments.');
end

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'ptx','en_nl','nl_char','nl_lim', 'en_bp','bp_ord','bp_att','bp_fcut','fs'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% no signal: don't waste time
if ~justgauge && (isempty(signal) || sum(signal)==0)
   if isempty(signal)
      critwarn('Empty input signal.');
   else
      critwarn('Input signal has zero energy.');
   end
   return
end 


% *******************************************************************************************************
% prepare bandpass (TX antenna and other filters)

% ... use zpk and cascaded filter to avoid instability due to coefficient rounding
if settings.en_bp
   [z,p,k] = cheby2(settings.bp_ord, settings.bp_att, settings.bp_fcut*2/settings.fs);
   [s,g] = zp2sos(z,p,k);
   Hd_bp = dfilt.df1tsos(s,g); 
   % freqz(Hd_bp, linspace(settings.bp_fcut(1), 2*settings.bp_fcut(2), 1000), settings.fs); title('TX BANDPASS');
else
   [s,g] = zp2sos([], [], 1);
   Hd_bp = dfilt.df1tsos(s,g);
end

% just return (overall) filter object?
if justgauge
   % add amplifier gain
   % ATTENTION: *) transmit power is set AFTER the filter
   %            *) normalization of input variance (additional gain) cannot be taken into account here
   [s_ps,g_ps] = zp2sos([], [], sqrt(settings.ptx / sum(abs(impz(Hd_bp)).^2)) );
   signal = dfilt.cascade(Hd_bp, dfilt.df1tsos(s_ps,g_ps));
   return
end


% *******************************************************************************************************
% power amplifier

% nonlinearity (if enabled)
if settings.en_nl
   % load characteristic
   pwramp = loadmat(settings.nl_char);
   % this works only if the characteristic is symmetrical => check (assume full symmetry; do not interp)
   if abs(mean(pwramp.nl_x)) > internalsettings.pwramp_symtol || abs(mean(pwramp.nl_y)) > internalsettings.pwramp_symtol
      criterr('Characteristic of power-amplifier has to be symmetrical around zero.');
   end
   % determine input range of signal and scale characteristic accordingly (char. is x=-1...1 and |y|<=1)
   signal_min = min(signal);
   signal_max = max(signal);
   pwramp.nl_x = (pwramp.nl_x/settings.nl_lim+1)/2 * (signal_max-signal_min) + signal_min;
   pwramp.nl_y = (pwramp.nl_y                +1)/2 * (signal_max-signal_min) + signal_min;   
   % warp (extrap. for borders)
   signal = interp1(pwramp.nl_x, pwramp.nl_y, signal, 'linear', 'extrap'); 
end


% *******************************************************************************************************
% filter, normalize and set power (also: compensate for power lost in filters)

% input bandpass (if enabled)
if settings.en_bp
   signal = filter(Hd_bp, signal);
end

% set power
signal = signal * sqrt(settings.ptx/var(signal)); 
