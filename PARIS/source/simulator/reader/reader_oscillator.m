% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% reader - oscillator
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
% version = reader_oscillator()
%    Just returns the version number (string).
% carrier = reader_oscillator(settings)
%    Creates a cosine, sine or complex exponential wave with variance = 1 (for exp: var(re)=var(im)=1 !)
%    and nonconstant frequency according to SETTINGS. Note that only the average variance is exactly 
%    equal to one (not normalized for each vector).
%
%
% ***** Interface definition *****
% function carrier = reader_oscillator(settings)
%    settings   struct containing all necessary parameters
%       .fcenter        center frequency in Hz
%       .fstddev        frequency standard deviation in Hz (Gaussian) ... phase noise
%       .astddev        amplitude standard deviation (center = 1; Gaussian)
%       .snr            SNR in dBc (to carrier, i.e., to 0dB)
%       .length         length of carrier wave signal in s
%       .fs             sampling frequency in Hz
%       .mode           (optional) {'sin', 'cos', 'exp'} create sine/cosine/exponential (default: cos)
%      
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

function carrier = reader_oscillator(settings)
version = 'beta 3.0';


% *******************************************************************************************************
% internal settings

% minimum oversampling rates (fs / (2*fcenter)) before ...
internalsettings.os_warn     = 8; % ... warning is issued
internalsettings.os_critwarn = 2; % ... critical warning is issued
internalsettings.os_err      = 1; % ... error (subsampling)


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   carrier = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'fcenter', 'fstddev', 'astddev', 'snr', 'length', 'fs'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% warning in case of low oversampling rate, error in case of subsampling
settings.os = settings.fs / (2*settings.fcenter);
if settings.os < internalsettings.os_err
   err('Oversampling rate (%.3f) is too low -> subsampling.', settings.os)
elseif settings.os < internalsettings.os_critwarn
   critwarn('Oversampling rate (%.3f) is close to critical.', settings.os);
elseif settings.os < internalsettings.os_warn
   warn('Low oversampling rate (%.3f); carrier peak amplitudes might vary significantly.', settings.os);
end

% optional parameter settings.mode
if ~isfield(settings, 'mode')
   settings.mode = 'cos';
end


% *******************************************************************************************************
% create carrier wave cos/sin/exp(w*t)
%  ... doing everything in one step would be more efficient, at the cost of (all) readability

% vector length (samples)
settings.length_s = round(settings.length * settings.fs);

% step1: time / frequency
if settings.fstddev == 0
   carrier = settings.fcenter/settings.fs * [0:1:settings.length_s-1]';
else
   carrier = settings.fcenter/settings.fs * [0:1:settings.length_s-1]' + ...
      + randn(settings.length_s, 1) * settings.fstddev/settings.fcenter;
end

% step2: function and amplitude instability
if settings.astddev == 0
   switch lower(settings.mode)
      case 'sin'
         carrier = sqrt(2) * sin(2*pi*carrier);
      case 'cos'
         carrier = sqrt(2) * cos(2*pi*carrier);
      case 'exp'
         carrier = sqrt(2) * exp(complex(0, 2*pi*carrier));
      otherwise
         err('Unsupported settings.mode = "%s".', settings.mode);
   end
else
   switch lower(settings.mode)
      case 'sin'
         carrier = sqrt(2) * (1 + (randn(settings.length_s, 1) * settings.astddev)) .* sin(2*pi*carrier);
      case 'cos'
         carrier = sqrt(2) * (1 + (randn(settings.length_s, 1) * settings.astddev)) .* cos(2*pi*carrier);
      case 'exp'
         carrier = sqrt(2) * (1 + (randn(settings.length_s, 1) * settings.astddev)) .* exp(complex(0, 2*pi*carrier));
      otherwise
         err('Unsupported settings.mode = "%s".', settings.mode);
   end
end

% step3: additive white gaussian noise
varnoise = 10^(-settings.snr/10);
if varnoise ~= 0
   if strcmpi(settings.mode, 'exp')
      carrier = carrier + complex(randn(settings.length_s, 1),randn(settings.length_s, 1))*sqrt(varnoise);
   else
      carrier = carrier + randn(settings.length_s, 1)*sqrt(varnoise);
   end
end

