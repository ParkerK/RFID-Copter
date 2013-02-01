% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% reader - receiver (input) stage (bandpass, power splitter)
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
% version = reader_receiver()
%    Just returns the version number (string).
% Hd = reader_receiver(settings)
%    Returns the overall (linear) filter object of this function (dfilt); valid for I and Q channel.
% [sig_i, sig_q] = reader_receiver(signal, settings)
%    Applies the reader input stage to SIGNAL and returns I and Q channel signals (after power splitter)
%    in a complex vector. Filtering is done by a cascaded second order IIR structure.
%
%
% ***** Interface definition *****
% function signal = reader_receiver(signal, settings)
%    signal        (optional) input signal
%    settings      struct containing settings
%       .en_bp        enable bandpass true/false
%       .iirord       input bandpass filter order
%       .att          stopband attenuation of input bandpass in dB
%       .fcut         input bandpass stopband edge frequencies in Hz
%       .fs           sampling frequency in Hz
%
%    signal  complex received signal (I+jQ)
%            -> alternatively: filter object (overall dfilt filter structure of this function) 
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
% - make filter a multi/narrowband bandpass (mfcw carrier pairs w. subsampling)
%
% *******************************************************************************************************

function signal = reader_receiver(varargin)
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
expected.reqfields = {'en_bp', 'iirord', 'att', 'fcut', 'fs'};
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
% generate receiver structure

% input bandpass
% ... use zpk and cascaded filter to avoid instability due to coefficient rounding
if settings.en_bp
   [z,p,k] = cheby2(settings.iirord, settings.att, settings.fcut*2/settings.fs);
   [s,g] = zp2sos(z,p,k);
   Hd_bp = dfilt.df1tsos(s,g); % freqz(Hd_bp, linspace(settings.fcut(1), 5*settings.fcut(2), 1000), settings.fs); title('RX BANDPASS');
else
   [s,g] = zp2sos([], [], 1);
   Hd_bp = dfilt.df1tsos(s,g);
end

% just return (overall) filter object?
if justgauge
   [s_ps,g_ps] = zp2sos([], [], 1/2); % power splitter
   signal = dfilt.cascade(Hd_bp, dfilt.df1tsos(s_ps,g_ps));
   return
end


% *******************************************************************************************************
% apply receiver structure

% input bandpass (if enabled)
if settings.en_bp
   signal = filter(Hd_bp, signal);
end

% power splitter / switch to complex
signal = complex(1/2*signal, 1/2*signal);
