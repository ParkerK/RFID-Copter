% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% reader - estimates the reader "analog" path (transmitter, receiver, demodulator) spectrum  
%          at certain frequencies
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
% version = reader_analogpathest()
%    Just returns the version number (string).
% hfreq = reader_analogpathest(freq, settings, txsettings, rxsettings, demodsettings)
%    Returns the frequency response of all reader analog path functions at frequencies FREQ (baseband!). 
%    Systematic and random deviations are also modeled (controlled by settings). Note that the errors are
%    added to the responses of all functions and the overall path independently, i.e. the overall 
%    frequency response is NOT the concatenation of all blocks if errors are added.
%
%
% ***** Interface definition *****
% hfreq = reader_analogpathest(freq, settings, rxsettings, demodsettings)
%    freq            baseband frequency vector in Hz (can be scalar)
%    settings        struct containing settings
%       .f0             carrier center frequency in Hz (the carrier that is used for demodulation)
%       .fs             sampling frequency in Hz
%       .err         error in estimates: ~Gaussian [systematic gain, random gain] ... no error: [1,0]
%                    [mult. factor (syst. err.), std. dev. factor (rand err.)]
%    txsettings      settings for READER_TRANSMITTER (see reader_transmitter.m for further details)
%    rxsettings      settings for READER_RECEIVER (see reader_receiver.m for further details)
%    demodsettings   settings for READER_DEMODULATOR (see reader_demodulator.m for further details)
%
%    hfreq   struct with frequency responses of all reader analog path functions at frequencies FREQU
%       .transmitter    response of READER_TRANSMITTER at frequency SETTINGS.F0 + FREQU
%       .receiver       response of READER_RECEIVER at frequency SETTINGS.F0 + FREQU
%       .demodulation   response of READER_DEMODULATOR at frequency FREQU
%       .all            combination of .receiver and .demodulation (overall response)
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

function hfreq = reader_analogpathest(freq, settings, txsettings, rxsettings, demodsettings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   hfreq = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'f0', 'fs', 'err'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% frequencies of interest
freq = [0, freq(:)']; % make sure freq is a row /vector/ (freqz also creates row vectors)


% *******************************************************************************************************
% get linear models of input stage and response @ freq

% transmitter: around carrier f0
hd.transmitter       = reader_transmitter(txsettings);
hfreq.transmitter    = freqz(hd.transmitter, settings.f0 + freq, settings.fs);
hfreq.transmitter(1) = [];

% receiver: around carrier f0
hd.receiver       = reader_receiver(rxsettings);
hfreq.receiver    = freqz(hd.receiver, settings.f0 + freq, settings.fs);
hfreq.receiver(1) = [];

% demodulator: around DC, but still @ fs
hd.demodulation       = reader_demodulation(demodsettings);
hfreq.demodulation    = freqz(hd.demodulation, freq, settings.fs);
hfreq.demodulation(1) = [];

% combine
hfreq.all = hfreq.transmitter .* hfreq.receiver .* hfreq.demodulation;


% *******************************************************************************************************
% introduce error ("estimate")

hfreq.transmitter  = ( settings.err(1) + randn(1,length(freq)-1)*settings.err(2) ) .* hfreq.transmitter;
hfreq.receiver     = ( settings.err(1) + randn(1,length(freq)-1)*settings.err(2) ) .* hfreq.receiver;
hfreq.demodulation = ( settings.err(1) + randn(1,length(freq)-1)*settings.err(2) ) .* hfreq.demodulation;
hfreq.all          = ( settings.err(1) + randn(1,length(freq)-1)*settings.err(2) ) .* hfreq.all;
