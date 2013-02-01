% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - obtain (tag) modulation signal timing
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
% timings = linktiming_tag()
%    Just returns the version number (string).
% timings = linktiming_tag(signal)
%    Returns a struct containing measured timings and an envelope analysis in samples. Also does some
%    protocol analysis (encoding, pilot tone). This function assumes a time-invariant envelope.
% timings = linktiming_tag(signal, settings)
%    Same functionality as "linktiming_tag(signal)", but with some user-controllable parameters. 
%
%
% ***** Function definition *****
% timings = linktiming_reader(signal)
%    signal     modulated carrier
%    settings   (optional) struct containing settings
%       .levels    [low, high] level thresholds (e.g. [0.1, 0.9] means 10%, 90%) ... values close to 0.5
%                  make the result more robust against overshoots at the cost of accuracy
%
%    info       struct containing the measured (estimated) values in samples
%       .am1           average high (modulated) amplitude
%       .am0           average low (unmodulated) amplitude
%       .os_r          [avg, std. dev.] max. amplitude of overshoots at rising edges (start of "high")
%       .os_f          [avg, std. dev.] max. amplitude of overshoots at falling edges (end of "high")
%       .us_r          [avg, std. dev.] min. amplitude of undershoots at rising edges (end of "low")
%       .us_f          [avg, std. dev.] min. amplitude of undershoots at falling edges (start of "low")
%       .tsub_s        [avg, std. dev., N] subcarrier HALFperiod in samples, N: sum of occurrence
%       .tinv_s        [avg, std. dev., N] subc. HALFperiod at inversion in samples, N: sum of occurrence
%       .tsub_s        [avg, std. dev., N] subc. HALFperiod at violation in samples, N: sum of occurrence
%       .tlf_s         subcarrier period in samples (refined 2*tsub_s)
%       .enc_m         encoding: M (1:FM0, 2,4,8:Miller)
%       .enc_cert      encoding: certainty of preamble correlation, i.e. ratio matching/nonmatching
%       .enc_trext     encoding: pilot tone detected (1/0); NaN if an invalid pilot signal was detected
%       .tf_s          [avg, std deviation] of fall times tf in samples
%       .tr_s          [avg, std deviation] of rise times tr in samples
%       .t0_s          length of unmodulated parts at [beginning, end] in samples (to 50% threshold)
%       .trans_s       position of transients (subbit borders) in signal
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
% ? make robust against noise/no signal
%   
% *******************************************************************************************************

function info = linktiming_tag(varargin)
version = 'beta 3.0';

% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   info = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% thresholds 
internalsettings.dc = eps; % signal is DC if (signal - mean) / mean < .dc 
% ... the closer to 0.5, the more robust against overshoots, at the cost of accuracy
internalsettings.levels = [0.1, 0.9]; % [low, high] thresholds

% demodulation
internalsettings.demod.order     =     4; % filter order (/2 because filtfilt is used)
internalsettings.demod.fcut      = 1/100; % normalized cutoff frequency (1: fs/2)
internalsettings.demod.att       =    40; % dB stopband attenuation  (/2 because filtfilt is used)
internalsettings.demod.transient =   250; % samples "transient" area at start/end of signal: removed
internalsettings.downsampling    =    50; % downsampling factor (tradeoff speed <-> resolution)

% preamble: FM0
internalsettings.pre0 = [1,1,-1,1,-1,-1,1,-1,-1,-1,1,1]'; % FM0, subcarrier level
% preamble: Miller (create using tag_encoding without data, don't add CRC)
%     deactivate warnings ("input vector of length zero")
global globalsettings;
globalsettings_old = globalsettings;
globalsettings.logging.warnings = 0;
%     get preambles
internalsettings.encoding.trext = 0; % no pilot tone
internalsettings.encoding.m = 2; % Miller M=2
internalsettings.pre2 = tag_encoding([], internalsettings.encoding, false);
internalsettings.encoding.m = 4; % Miller M=4
internalsettings.pre4 = tag_encoding([], internalsettings.encoding, false);
internalsettings.encoding.m = 8; % Miller M=8
internalsettings.pre8 = tag_encoding([], internalsettings.encoding, false);
%     reset warning setting
globalsettings = globalsettings_old;
%     remove end-of-signaling 
internalsettings.pre2 = internalsettings.pre2(1:end- 4);
internalsettings.pre4 = internalsettings.pre4(1:end- 8);
internalsettings.pre8 = internalsettings.pre8(1:end-16);


% *******************************************************************************************************
% input parameter checks / prepare input parameters

% number of input variables
if nargin < 1
    criterr('Not enough input arguments.')
elseif nargin == 1
   signal          = varargin{1};
   settings.levels = internalsettings.levels;
elseif nargin == 2
   signal   = varargin{1};
   settings = varargin{2};
else
    criterr('Too many input parameters.')  
end

% check contents of settings
%     prepare required data
expected.name = 'settings';
expected.reqfields = {'levels'};
%     check
errortext = contentcheck(settings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% length of input vector
if isempty(signal)
    criterr('Length of input signal is zero.')
end

% just DC
if all(signal - mean(signal) > internalsettings.dc*mean(signal))
   critwarn('Constant input signal. Returning an empty struct.');
   info = struct();
   return
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
%     ... the demodulated signal always starts with a rising edge

% get modulation amplitudes A,B (first estimate)
am_1 = max(signal); % A
am_0 = min(signal); % B

% get intersection with low and high threshold (not 50% because of better performance at start/end)
am_h = am_0 + (am_1 - am_0) * settings.levels(2);
am_l = am_0 + (am_1 - am_0) * settings.levels(1);
int_h = findzeros(signal - am_h);
int_l = findzeros(signal - am_l);
 % add "start of window" and "end of window"
int_h = [1; int_h; length(signal)];
int_l = [1; int_l; length(signal)];

% now: int%%(even) -> int%%(odd) : high "1"
%      int%%(odd)  -> int%%(even): low "0"
% ... take index 1/2 between crossings and recalculate am_1 and am_0 values
%     (median because border regions might be cut off for small t0)
am_1 = median(signal( round( (int_h(2:2:end-1) + int_h(3:2:end)) / 2 ) )); % attention: start/end included
am_0 = median(signal( round( (int_l(1:2:end)   + int_l(2:2:end)) / 2 ) )); % attention: start/end included

% calculate low, 50%, high values
am_h = am_0 + (am_1 - am_0) * settings.levels(2);
am50 = am_0 + (am_1 - am_0) * 0.5;
am_l = am_0 + (am_1 - am_0) * settings.levels(1);

% get intersections with thresholds
%     int50 could also be done with filtfilt(b,a, abs(diff(signal)));
int_h = findzeros(signal - am_h);
int50 = findzeros(signal - am50); 
int_l = findzeros(signal - am_l);


% *******************************************************************************************************
% amplitude calculations I (not prone to overshoot problems)

% high/low amplitudes
info.am1 = am_1;
info.am0 = am_0;


% *******************************************************************************************************
% timing calculations (everything in samples @ fs)
%     ... the demodulated signal always starts with a rising edge

% unmodulated parts at start/end (coarse estimates); compensate for downsampling
% ...  also correct for linear slope: low-high => 0-100 (estimated trf are only 80%)
info.t0_s = [int_l(1) - 1 - (int_h(1)-int_l(1)),...
             length(signal) - int_l(end) - (int_l(end)-int_h(end))] *...
             internalsettings.downsampling;

% subcarrier transients
dur = diff(int50);
          
% add transient positions to info struct, correct for downsampling
info.tran_s = int50 * internalsettings.downsampling;

% now there should be 2 or 3 classes: t_subcarrier, inversion = 2*t_sub, optional: violation = 3*t_sub
% ... use a histogram to separate these classes ...
[n, t_s] = hist(dur, min(dur):1:max(dur));
n   = sparse(n);
% ... and calculate expected values (mean, variance); compensate for downsampling
%    ... if there is anything in between 1/3 and 3/4 (3/4 because violations are much less likely than
%        subcarrier period => less jitter problems): 3 subclasses
if full(sum(n(round(end/3):round(3/4*end)))) ~= 0 
   info.tsub_s = [[1/sum(n(1:round(end/3))) * sum(t_s(1:round(end/3)) .* n(1:round(end/3))),... % E{n}
      sqrt(var(t_s(1:round(end/3)), n(1:round(end/3))))] * internalsettings.downsampling,... % std. deviation
      sum(n(1:round(end/3)))]; % sum of occurrence
   info.tinv_s = [[1/sum(n(round(end/3)+1:round(end*3/4)))*sum(t_s(round(end/3)+1:round(end*3/4)).*n(round(end/3)+1:round(end*3/4))),...
      sqrt(var(t_s(round(end/3)+1:round(end*3/4)), n(round(end/3)+1:round(end*3/4))))] * internalsettings.downsampling,...
      sum(n(round(end/3)+1:round(end*3/4)))];
   info.tvio_s = [[1/sum(n(round(end*3/4)+1:end))*sum(t_s(round(end*3/4)+1:end).*n(round(end*3/4)+1:end)),...
      sqrt(var(t_s(round(end*3/4)+1:end), n(round(end*3/4)+1:end)))] * internalsettings.downsampling,...
      sum(n(round(end*3/4)+1:end))];
else % if not: 2 classes
   info.tsub_s = [[1/sum(n(1:round(end/2))) * sum(t_s(1:round(end/2)) .* n(1:round(end/2))),... % E{n}
      sqrt(var(t_s(1:round(end/2)), n(1:round(end/2))))] * internalsettings.downsampling,... % std. deviation
      sum(n(1:round(end/2)))];  % sum of occurrence
   info.tinv_s = [[1/sum(n(round(end/2)+1:end)) * sum(t_s(round(end/2)+1:end) .* n(round(end/2)+1:end)),...
      sqrt(var(t_s(round(end/2)+1:end), n(round(end/2)+1:end)))] * internalsettings.downsampling,...
      sum(n(round(end/2)+1:end))];
   info.tvio_s = [NaN, NaN, 0];
end
%     check if the results make sense (in case of protocol violations: only one or more than 3 classes)
if    info.tinv_s(1) < 1.5*info.tsub_s(1) || info.tinv_s(1) > 2.5*info.tsub_s(1) ||...
      info.tvio_s(1) < 2.5*info.tsub_s(1) || info.tvio_s(1) > 3.5*info.tsub_s(1)
   critwarn('Unable to detect timing; possible protocol violation. Assuming only subcarrier is present.');
   info.tsub_s = [[1/sum(n)*sum(t_s.*n), sqrt(var(t_s, n))]*internalsettings.downsampling, sum(n)];
   info.tvio_s = [NaN, NaN, 0];
   info.tinv_s = [NaN, NaN, 0];
end
%     try calculate LF period in samples (refined 2*tsub_s)
info.tlf_s = (int50(end)-int50(1))*internalsettings.downsampling * 2 /... % overall time * 2 (=> period!)
             (info.tsub_s(3) + 2*info.tinv_s(3) + 3*info.tvio_s(3)); % divided by number of subsections


% check if level intersections make sense
if length(int_h) ~= length(int50) || length(int50) ~= length(int_l)
   critwarn('Massive overshoots detected. Timing/amplitude/protocol detection will fail; skipping.');
   info.massive_os = true;
   % set dummy values
   info.tf_s = [NaN, NaN];
   info.tr_s = [NaN, NaN];
   info.enc_m = NaN;
   info.enc_cert = NaN;
   info.enc_trext = NaN;
   info.os_r = [NaN, NaN];
   info.os_f = [NaN, NaN];
   info.us_r = [NaN, NaN];
   info.us_f = [NaN, NaN];
   % and skip all following calculations
   return
else
   info.massive_os = false;
end

% rising / falling edges; compensate for downsampling
% ... also correct for linear slope: low-high => 0-100
tr = int_h(1:2:end) - int_l(1:2:end);
tf = int_l(2:2:end) - int_h(2:2:end);
info.tf_s = [mean(tf), std(tf)] * internalsettings.downsampling/(settings.levels(2)-settings.levels(1));
info.tr_s = [mean(tr), std(tr)] * internalsettings.downsampling/(settings.levels(2)-settings.levels(1));


% *******************************************************************************************************
% amplitude calculations II (prone to overshoot problems)
%     ... the demodulated signal always starts with a rising edge

% envelope: overshoots/undershoots (does not have to be very performant)
%     overshoots (@ modulated)
os_r = zeros(length(int_h)/2, 1); % @ rising edge
os_f = zeros(length(int_h)/2, 1); % @ falling edge
for i = 1 : length(int_h)/2
   os_r(i) = max(signal( int_h(2*i-1):round((int_h(2*i-1)+int_h(2*i))/2)   ));
   os_f(i) = max(signal( round((int_h(2*i-1)+int_h(2*i))/2):int_h(2*i) ));
end
info.os_r = [mean(os_r), std(os_r)];
info.os_f = [mean(os_f), std(os_f)];
%     undershoots (@ unmodulated)
int_l = [1; int_l; length(signal)]; % also include first/last transient
us_r = zeros(length(int_l)/2, 1); % @ rising edge
us_f = zeros(length(int_l)/2, 1); % @ falling edge
for i = 1 : length(int_l)/2
   us_r(i) = min(signal( int_l(2*i-1):round((int_l(2*i-1)+int_l(2*i))/2)   ));
   us_f(i) = min(signal( round((int_l(2*i-1)+int_l(2*i))/2):int_l(2*i) ));
end
int_l([1,end]) = []; % limit to transients again
info.us_r = [mean(us_r), std(us_r)];
info.us_f = [mean(us_f), std(us_f)];


% *******************************************************************************************************
% try do detect encoding

% quantize signal using the detected subcarrier halfperiod 
qsig = []; % {-1,1} (better suited for correlation
for i = 1 : length(int50) - 1
   dt_s = (int50(i+1) - int50(i)) * internalsettings.downsampling;
   sbit = (mean(signal(int50(i):int50(i+1))) > am50)*2-1;
   if dt_s < 1.5 * info.tlf_s(1)/2 % normal halfcarrier period
      qsig = [qsig; sbit];
   elseif dt_s > 1.5 * info.tlf_s(1)/2 && dt_s < 2.5 * info.tlf_s(1)/2 % inversion
      qsig = [qsig; sbit; sbit];
   elseif dt_s > 2.5 * info.tlf_s(1)/2 && dt_s < 3.5 * info.tlf_s(1)/2 % protocol violation
      qsig = [qsig; sbit*ones(3,1)];
   else % this should not happen
      critwarn('Unable to detect encoding; possible protocol violation.');
      qsig = [qsig; NaN];
   end
end

% correlate with preambles, select maxima
pre_corr = [max(xcorr(qsig, internalsettings.pre0))/length(internalsettings.pre0),...
            max(xcorr(qsig, internalsettings.pre2))/length(internalsettings.pre2),...
            max(xcorr(qsig, internalsettings.pre4))/length(internalsettings.pre4),...
            max(xcorr(qsig, internalsettings.pre8))/length(internalsettings.pre8)];
%     determine best match, set encoding (m) and preamble accordingly
[dummy, ind] = max(pre_corr);
switch ind
   case 1 % FM0
      info.enc_m = 1;
      preamble = internalsettings.pre0;
   case 2 % Miller, M=2
      info.enc_m = 2;
      preamble = internalsettings.pre2;
   case 3 % Miller, M=4
      info.enc_m = 4;
      preamble = internalsettings.pre4;
   case 4 % Miller, M=5
      info.enc_m = 8;
      preamble = internalsettings.pre8;
   otherwise
end
%     pilot tone?
[match, ind] = max(fftshift(xcorr(qsig, preamble)));
ind = mod(ind, length(qsig)-1);
info.enc_cert = match/length(preamble); % certainty (ratio of matching of preamble)
if ind == 1 % preamble at beginning => no pilot tone
   info.enc_trext = 0;
else
   if info.enc_m == 1 % FM0: 12 subcarrier periods
      if ind == 24 % these are halfperiods
         info.enc_trext = 1;
      else
         warn('Unable to determine if a pilot tone is present (FM0).');
         info.enc_trext = NaN;
      end
   else % Miller: 12 (16 - 4 already in preamble) carrier periods
      if ind == (16-4) * 2 * info.enc_m % these are halfperiods
         info.enc_trext = 1;
      else
         warn('Unable to determine if a pilot tone is present (Miller).');
         info.enc_trext = NaN;
      end
   end
end



return
% % *******************************************************************************************************
% % DEBUG
% close all
% 
% % signals and levels (Warning: Assumes EPCGlobal Protocol)
% figure; hold on;
% plot(signal, 'b-');
% plot(ones(size(signal))*am_1, 'g-');
% plot(ones(size(signal))*am_h, 'g-');
% plot(ones(size(signal))*am50, 'g--');
% plot(ones(size(signal))*am_l, 'g--');
% plot(ones(size(signal))*am_0, 'g--');
% plot(int_h, ones(size(int_h))*am_h, 'ko');
% plot(int50, ones(size(int50))*am50, 'ko');
% plot(int_l, ones(size(int_l))*am_l, 'ko');
% plot(info.tlf_s(1)/2/internalsettings.downsampling*[0.5:1:length(qsig)-0.5]+int50(1),...
%    am_0 + (qsig+1)/2*(am_1-am_0), 'ro'); % this will not work for very long signals
% plot(info.tlf_s(1)/2/internalsettings.downsampling*[0:1:length(qsig)]+int50(1),...
%    am50, 'rx'); % this will not work for very long signals
% hold off; grid on; axis tight;
% 
% % timing classification
% if ~isnan(info.tvio_s)
%    pdf1 = normpdf(t_s*internalsettings.downsampling, info.tsub_s(1), info.tsub_s(2));
%    pdf1 = pdf1/max(pdf1)*max(n(1:round(end/3)));
%    pdf2 = normpdf(t_s*internalsettings.downsampling, info.tinv_s(1), info.tinv_s(2));
%    pdf2 = pdf2/max(pdf2)*max(n(round(end/3)+1:round(end*3/4)));
%    pdf3 = normpdf(t_s*internalsettings.downsampling, info.tvio_s(1), info.tvio_s(2));
%    pdf3 = pdf3/max(pdf3)*max(n(round(end*3/4)+1:end));
% else
%    pdf1 = normpdf(t_s*internalsettings.downsampling, info.tsub_s(1), info.tsub_s(2));
%    pdf1 = pdf1/max(pdf1)*max(n(1:round(end/2)));
%    pdf2 = normpdf(t_s*internalsettings.downsampling, info.tinv_s(1), info.tinv_s(2));
%    pdf2 = pdf2/max(pdf2)*max(n(round(end/2)+1:end));
%    pdf3 = nan(size(pdf1));
% end
% figure; hold on;
% stem(t_s*internalsettings.downsampling, n);
% plot(t_s*internalsettings.downsampling, pdf1, 'r');
% plot(t_s*internalsettings.downsampling, pdf2, 'g');
% plot(t_s*internalsettings.downsampling, pdf3, 'k');
% stem((t_s(1)+[0, 1/4, 1/3, 1/2, 2/3, 3/4, 1]*(t_s(end)-t_s(1)))*internalsettings.downsampling, ones(7,1)*max(n), 'k--', 'Marker','none');
% hold off; xlabel('t [samples @ fs]');
% 
% % preamble correlation
% figure; plot(circshift(fftshift(xcorr(qsig, internalsettings.pre0))/length(internalsettings.pre0), 1))
% figure; plot(circshift(fftshift(xcorr(qsig, internalsettings.pre2))/length(internalsettings.pre2), 1))
% figure; plot(circshift(fftshift(xcorr(qsig, internalsettings.pre4))/length(internalsettings.pre4), 1))
% figure; plot(circshift(fftshift(xcorr(qsig, internalsettings.pre8))/length(internalsettings.pre8), 1))


% *******************************************************************************************************
% REPLACED PARTS

% info.t0_s = [int_l(1) - 1 - (int_h(1)-int_l(1))*...
%              settings.levels(1)/(settings.levels(2)-settings.levels(1)),...
%              length(signal) - int_l(end) - (int_l(end)-int_h(end))*...
%              settings.levels(1)/(settings.levels(2)-settings.levels(1))] *...
%              internalsettings.downsampling;
