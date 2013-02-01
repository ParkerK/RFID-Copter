% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - obtain (reader) modulation signal timing
%
% Note that this function is not robust againgst noise, degraded signals or protocol changes/problems. 
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
% timings = linktiming_reader()
%    Just returns the version number (string).
% timings = linktiming_reader(signal)
%    Returns a struct containing measured timings in samples. Assumes cosine rolloff and
%    time-invariant timings. Please note that very small times (esp. t0, tr, tf) are measured with high
%    variance. Timings close to zero (t0, tr, tf, ...) are not supported and may cause errors.
%
%
% ***** Function definition *****
% timings = linktiming_reader(signal)
%    signal   modulated carrier
%
%    timings   struct containing the measured timings in samples
%       .moddepth      [avg, std deviation] modulation depth in percent
%       .delimiter_s   length of delimiter (preamble/framesync)
%       .tari_s        Tari (preamble/framesync)
%       .rtcal_s       RTcal (preamble/framesync)
%       .trcal_s       TRcal (preamble, NaN for framesync)
%       .data          decoded data bits
%       .data0_s       [avg length, std deviation] of data-0 symbols
%       .data1_s       [avg length, std deviation] of data-1 symbols
%       .x_s           [avg length, std deviation] of x (additional length of data-1)
%       .pw_s          [avg length, std deviation] of pulse widths PW
%       .tf_s          [avg length, std deviation] of fall times tf (10-90%)
%       .tr_s          [avg length, std deviation] of rise times tr (10-90%)
%       .t0_s          length of unmodulated parts at beginning / end (90%)
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
% = implement ripple measurement
% = time-variant moddepth measurement
% - add monotony checks (avoid problems for small timings, e.g. trf=0)
% - improve 10/50/90% value generation (time variant?)
% - there might be some systematic error for t0_s
% ? stddev for t0?
%   
% *******************************************************************************************************

function timings = linktiming_reader(signal)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   timings = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings
internalsettings.demod.order     = 4;    % filter order (/2 because filtfilt is used)
internalsettings.demod.fcut      = 1/10; % normalized cutoff frequency (1: fs/2)
internalsettings.demod.att       = 40;   % dB stopband attenuation  (/2 because filtfilt is used)
internalsettings.demod.transient = 250;  % samples "transient" area at start/end of signal: removed
internalsettings.downsampling    = 25;   % downsampling factor (tradeoff speed <-> resolution)


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% number of input variables
if nargin < 1
    criterr('Not enough input arguments.')
end

% length of input vector
if isempty(signal)
    criterr('Length of input signal is zero.')
end


% *******************************************************************************************************
% straightforward demodulation and downsampling

% design filter (Chebycheff type 2: stopband ripple)
[b,a] = cheby2(internalsettings.demod.order, internalsettings.demod.att, internalsettings.demod.fcut);

% rectifier and filter => simple AM demod
signal = filtfilt(b,a, abs(signal));

% delete transients and downsampling
signal = signal( internalsettings.demod.transient:internalsettings.downsampling:end-internalsettings.demod.transient );


% *******************************************************************************************************
% amplitude and timing detection

% get modulation amplitudes A,B (first estimate)
am_1 = max(signal); % A
am_0 = min(signal); % B

% get intersection with 10 and 90% threshold (not 50% because of better performance at start/end)
am90 = am_0 + (am_1 - am_0) * 0.9;
am10 = am_0 + (am_1 - am_0) * 0.1;
int90 = findzeros(signal - am90);
int10 = findzeros(signal - am10);
 % add "start of window" and "end of window"
int90 = [1; int90; length(signal)];
int10 = [1; int10; length(signal)];

% now: int%%(odd)  -> int%%(even): high (A) ; 1->2 , 3->4 , 5->6 , ...
%      int%%(even) -> int%%(odd) : low (B)  ; 2->3 , 4->5 , ...
% ... take index 1/2 between crossings
% and recalculate am_1 and am_0 values (median because border regions might be cut off for small t0)
am_1 = median(signal( round( (int90(1:2:end)   + int90(2:2:end)) / 2 ) )); % A
am_0 = median(signal( round( (int10(2:2:end-1) + int10(3:2:end)) / 2 ) )); % B

% calculate 10, 50, 90% values
am90 = am_0 + (am_1 - am_0) * 0.9;
am50 = am_0 + (am_1 - am_0) * 0.5;
am10 = am_0 + (am_1 - am_0) * 0.1;

% get intersections with thresholds, correct for downsampling
%     int50 could also be done with filtfilt(b,a, abs(diff(signal)));
int90 = findzeros(signal - am90) * internalsettings.downsampling;
int50 = findzeros(signal - am50) * internalsettings.downsampling; 
int10 = findzeros(signal - am10) * internalsettings.downsampling;


% *******************************************************************************************************
% amplitude calculations

timings.moddepth = [(am_1-am_0) / am_1 * 100 , -1];


% *******************************************************************************************************
% timing calculations (everything in samples @ fs)
%     ... the demodulated signal always starts with a falling edge (modulated)

% length of delimiter
timings.delimiter_s = int50(2) - int50(1);

% bit boundaries (including preamble)
%     difference between rising edges = int%%(even)
bits = diff(int50(2:2:end));

% preamble / framesync
timings.tari_s  = bits(1);
timings.rtcal_s = bits(2);
if bits(3) > timings.rtcal_s
   timings.trcal_s = bits(3);
   datastartindex = 4;
else
   timings.trcal_s = NaN;
   datastartindex = 3;
end

% data symbol duration
data   = [];
data_0 = [];
data_1 = [];
for i = datastartindex : length(bits)
   data = [data, bits(i) > 1.25 * timings.tari_s];
   if ~data(i-datastartindex+1) % data-0 (tari)
      data_0 = [data_0; bits(i)];
   else % data-1 (tari+x = 1.5...2 tari)
      data_1 = [data_1; bits(i)];
   end
end
timings.data  = data;
timings.data0_s = [mean(data_0), std(data_0)];
timings.data1_s = [mean(data_1), std(data_1)];
timings.x_s = [timings.data1_s(1)-timings.data0_s(1),...
   sqrt(timings.data1_s(2)^2+timings.data0_s(2)^2)]; % uncorrelated

% pulse width PW
%     ... falling->rising edge starting after delimiter
pws = int50(4:2:end) - int50(3:2:end);
timings.pw_s = [mean(pws), std(pws)];

% rising / falling edges
%     falling: odd indices
tf = int10(1:2:end) - int90(1:2:end);
tr = int90(2:2:end) - int10(2:2:end);
timings.tf_s = [mean(tf), std(tf)];
timings.tr_s = [mean(tr), std(tr)];

% coarse estimation of t0 (at start and end); compensate downsampling
timings.t0_s = (length(signal) * internalsettings.downsampling -...
   (int90(end)-int90(1)) ) / 2; % (overall length - modulated) / 2

% compensate for cos rolloff: cos 90-10% measured but 100-0% modulated
%     rise/fall time
trf_factor = (acos(-0.8) - acos(0.8)) / pi;
timings.tf_s = timings.tf_s / trf_factor;
timings.tr_s = timings.tr_s / trf_factor;
%     t0 (remove 100-90% of first falling slope)
timings.t0_s = timings.t0_s - timings.tf_s(1) * acos(0.8)/pi;


return
% *******************************************************************************************************
% DEBUG

% ind_A = round((int90(1:2:end)+int90(2:2:end))/2);
% ind_B = round((int10(2:2:end-1)+int10(3:2:end))/2);
% 
% figure; hold on;
% plot(signal, 'b-');
% plot(ones(size(signal))*am_1, 'g');
% plot(ones(size(signal))*am90, 'g');
% plot(ones(size(signal))*am50, 'g');
% plot(ones(size(signal))*am10, 'g');
% plot(ones(size(signal))*am_0, 'g');
% plot(ind_A/internalsettings.downsampling, signal(round(ind_A/internalsettings.downsampling)), 'ro');
% plot(ind_B/internalsettings.downsampling, signal(round(ind_B/internalsettings.downsampling)), 'bo');
% plot(int90/internalsettings.downsampling, ones(size(int90))*am90, 'ko');
% plot(int50/internalsettings.downsampling, ones(size(int50))*am50, 'ko');
% plot(int10/internalsettings.downsampling, ones(size(int10))*am10, 'ko');
% hold off; grid on;
