% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - FIR filter for sparse impulse responses and long input vector lengths
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
% version = sparse_fir()
%    Just returns the version number (string).
% output = sparse_fir(delay_s, gain, input, settings)
%    Calculates the filtered INPUT for an impulse response GAIN(DELAY_S), where GAIN is the magnitude of
%    the impulse response at discrete delay DELAY_S. The function checks the density (ratio of nonzero
%    elements) in GAIN an selects the most performant implementation based on a characteristic file. 
%    Note that the overhead for determining the fastest method can be considerable and even exceed 
%    the runtime of Matlab's filter() for short input signal lengths (< 1e3 samples).  Use filter() for
%    short input signals.
%    
%
%
% ***** Interface definition *****
% function output = sparse_fir(delay_s, gain, input, settings)
%    delays_s       discrete delays of filter impulse response
%    gain           magnitude of filter impulse response at delays delay_s 
%    input          input signal (to be filtered)
%    settings       struct containing settings
%       .charfile      filename for performance characteristic (tipping point)
%
%    output   filtered input signal
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
% = too slow (too much overhead)
% ? overhead check (very short input signals/impulse responses)
%
% *******************************************************************************************************

function output = sparse_fir(delay_s, gain, input, settings)

% % *******************************************************************************************************
% % input parameter checks / prepare input parameters
% 
% % check contents of settings
% %     prepare required data
% expected.name = 'settings';
% expected.reqfields = {'charfile'};
% %     check
% errortext = contentcheck(settings, expected);
% %     output
% if ~isempty(errortext)
%    err('Incomplete settings\n%s', errortext);
% end


% *******************************************************************************************************
% load performance characteristic, determine tipping point for filter methods

% load characteristic
pchar = loadmat(settings.charfile);

% find tipping point for density (interpolation in log scale)
%     input signal/impulse response length
len_i = length(input);
len_f = length(gain);
%     input signal length out of range
if len_i < pchar.len_i(1) && len_i > pchar.len_i(end)
   ind_f      = interp1(pchar.len_f_log10, pchar.sxing,          log10(len_f), 'nearest');
   tipping_pt = interp1(pchar.len_i_log10, pchar.sxing(:,ind_f), log10(len_i), 'linear', 'extrap');
%     impulse response length out of range
elseif len_f < pchar.len_f(1) && len_f > pchar.len_f(end)
   ind_i      = interp1(pchar.len_i_log10, pchar.sxing,          log10(len_i), 'nearest');
   tipping_pt = interp1(pchar.len_f_log10, pchar.sxing(ind_i,:), log10(len_f), 'linear', 'extrap');
%     inside ranges
else 
   [x, y] = meshgrid(pchar.len_f_log10, pchar.len_i_log10);
   tipping_pt = interp2(x,y, pchar.sxing, log10(len_i),log10(len_f), 'linear', NaN);
end


% *******************************************************************************************************
% filter

% loop faster (low density)
if nnz(delay_s)/length(delay_s) > tipping_pt
   output = zeros(size(input));
   for k = 1 : length(delay_s)
      output(1+delay_s(k):len_i) = output(1+delay_s(k):len_i) + input(1:end-delay_s(k)) * gain(k);
   end

% filter faster (high enough density); reconstruct full gain
else
   gain_full = zeros(1+delay_s(end), 1); % first term in numerator of FIR corresponds to delay_s=0
   gain_full(1+delay_s) = gain;
   output = filter(gain_full, 1, input);
end
