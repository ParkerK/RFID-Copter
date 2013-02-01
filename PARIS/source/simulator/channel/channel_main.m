% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% channel - main (time-invariant, XIXO)
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
% version = channel_main()
%    Just returns the version number (string).
% [output, statistics] = channel_main(input, settings)
%    Applies all channel models, adds noise at the receivers, and returns the result in OUTPUT.
%    Smallscale statistics and channel impulse reponses are returned in STATISTICS. This 
%    function supports and arbitrary number of inputs and outputs (XIXO).
%    The vectors in INPUT can be of arbitrary length; the length of all vectors in OUTPUT is identical
%    and equal to the maximum vector length in INPUT. The function assumes that noise settings depend
%    only on the receiver (identical across all transmitters for one receiver).
%    If the function is set to probe mode (settings.probemode = true) all outputs are empty vectors.
%
%    Note that the function sets globalsettings.core.mat_stubborn temporarily to true for loading of
%    sparse FIR performance characteristic.
% 
%    See also CHANNEL_LARGE, CHANNEL_DIRECTIVITY, CHANNEL_SURFACE, CHANNEL_SMALL, and CHANNEL_NOISE.
%
%
% ***** Global Variables *****
% globalsettings
%    .core.sfir_char       filename for sparse FIR performance characteristic
%    .core.mat_stubborn    ignore missing files; return empty struct if file is missing
%
%
% ***** Interface definition *****
% function output = channel_main(input, settings)
%    input       input vector cell array (one cell per transmitter)
%    settings    cell array of structs containing all necessary parameters (#TX times #RX)
%       .probemode       only record channels if true; do not calculate output
%       .type            channel type {'room': directivity only for LOS, 'outdoor': directivity for all paths)
%       .virtual         (optional) struct that defines this to be a virtual connection (mirrored TX/RX)
%          .vtx             true if the transmitter is virtual, false otherwise
%          .vrx             true if the receiver is virtual, false otherwise
%          .otx             originating transmitter for this virtual transmitter
%          .orx             originating receiver for this virtual receiver
%          .dim             dimension (normal) of the surface that created a virtual
%          .refl            number of reflections that lead to this virtual transmitter
%          .surface         (optional) link to the surface that created the virtual connection
%          .invsurf         (optional) switch to invert the linked surface ("everywhere except here")
%       .ucg             (optional) user controlled gain factor (can be used for attenuation due to multiple reflections)
%       .direct          (optional) struct: predefined large-scale gain/delay
%          .gain            gain factor
%          .delay_s         delay in samples
%       .largescale      struct containing large-scale model parameters
%          .type            type of large-scale model {'log-dist'}
%          .dist            distance/distances in m (supports up to two dimensions)
%          .pl              path-loss factor
%          .c               speed of light in m/s
%          .f0              "center" frequency/frequencies for attenuation in Hz
%          .fs              sampling frequency in Hz
%       .directivity     struct containing antenna directivity model parameters
%          .txant           filename of TX antenna characteristic (isotropic if empty)
%          .txrot           rotation of TX antenna [azimuth, elevation] ([0,0]: maximum in positive x-axis)
%          .rxant           filename of RX antenna characteristic (isotropic if empty)
%          .rxrot           rotation of RX antenna [azimuth, elevation] ([0,0]: maximum in positive x-axis)
%          .dir             direction of channel [azimuth, elevation] (from TX antenna)
%       .surfaces       struct containing substructs with definitions of one reflection each
%          .example        one such substruct, for example (e.g., definition of a wall)
%             .mode           mode for the surface: 0="none/off", 1="transmit", 2="reflect"
%             .vtx            this surface is covered by virtual transmitters
%             .blur           factor 0...1 blurring of borders (0: none, 1: full)
%                             also see INTERNALSETTINGS within this function
%             .n2             refractive index of reflective surface (typ n2>1)
%             .dim            dimension of the surface's normal vector (e.g. 3 means z in [x,y,z]) 
%             .pol_dim        dimension of polarization vector (electric field)  (e.g. 3 means z in [x,y,z]) 
%             .aoi            angle of incidence (0<=aoi<=90) in degree
%             .poi2e          dist. of the intersection w.t. surf. plane to the closest edge (<0: inside)
%             .largescale     large-scale model parameters for this reflection; see settings.largescale
%             .directivity    antenna directivity model parameters for this reflection; see settings.directivity
%       .smallscale      struct containing small-scale model parameters
%          .on              switch small-scale model on/off
%          .det             return the average PDP (true) or random NLOS paths with the same average PDP (false)
%          .seed            random seed for non-deterministic setup (set to NaN for normal operation)
%          .ensembles       number of independent ensembles for settings.det=false
%          .maxrays         maximum number of paths (not including the LOS path, i.e. maxrays=0 means LOS only)
%          .maxiter         maximum iterations for getting the rms delay spread right
%          .k               Ricean K-factor in dB
%          .trms            RMS delay spread in s
%          .bw              one-sided bandwidth for channel (with oversampling) in Hz
%          .fres            initial frequency resolution in Hz (might be increased to meet .trms)
%          .eps_k           relative tolerable error in K-factor (linear)
%          .eps_trms        relative tolerable error in RMS delay spread
%          .fs              sampling frequency in Hz
%       .noise           struct containing noise parameters
%          .type            type of noise {'wgn' white Gaussian, otherwise: no noise}
%                           (use 'off' to skip the call of CHANNEL_NOISE ... faster)
%          .n0              noise density [dBm, per ? Hz], e.g., [-50, 1000] = -50 dBm/kHz
%          .frxs            eceiver sampling frequency in Hz
%          .fs              sampling frequency in Hz
%       .c               speed of light in m/s
%       .fs              sampling frequency in Hz
%
%    output            output vector cell array (one cell per receiver)
%    statistics        cell array (#TX, #RX) ofs structs containing statistics returned by CHANNEL_SMALL 
%                      (if smallscale model is on), plus
%       .pdp_los          LOS component of power-delay-profile
%       .pdp_sum          sum of power-delay-profile
%       .k                K-factor in dB
%       .trms             RMS delay spread in s
%       .tmin             minimum delay in s
%       .tmax             maximum delay in s
%       .cir_d            channel impulse response: delay vector in s
%       .cir_g            channel impulse response: gain vector
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
% ! directivity gain for small-scale part in "outdoor" mode should only contain the TX gain
% ! spots with extremely high k-factor in results: NLOS parts removed for some reason?
% = revise parameter checks / filter part, ntx/nrx structure, txnum/rxnum workaround, CIR merge operation
% = revise surface / virtual transmitter handling
%     + use full raytracing to calculate reflection coefficients of VTX
%     + link VTX to all surfaces used in the creation; only switch off those
%     ? absorbtion of surfaces
%     ? right now refl. are not blocked by other surfaces
% - ultrawideband largescale channel
% ? mute version_system for large amount of channels (get versions only once)
% ? making entries in settings inactive by additional field in settings might be faster than deleting the
%   corresponding entry
% *******************************************************************************************************


function [output, statistics] = channel_main(input, settings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   output = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% warnings and messages
internalsettings.maxpaths_warn    = 250; % maximum number of (overall) paths before a warning is issued

% sparsity before loop implementation is selected (in case n performance characteristic is available)
internalsettings.default_xing_pt = 0.98; % <=0: always loop, >=1: always filter

% modes for surfaces (have to be identical to settings in CHANNEL_SURFACE and CHANNEL_NEWPOS)
internalsettings.surf_off  = 0; % neither reflection, nor transmission => surface inactive
internalsettings.surf_tran = 1; % "reflect" mode
internalsettings.surf_refl = 2; % "transmit" mode

% exponential rolloff close to surface edge (do not switch surface, but remove smoothly => diffraction)
internalsettings.surf_mingain =  1e-6; % minimum gain factor before gain == 0 and surface set to "off" (performance)
%     reflection mode (up round edges)
internalsettings.surf_refl_exp     =  8; % times wavelength @ fs transient zone width (cf. radar cross sections)
internalsettings.surf_refl_shift   =  1; % times wavelength @ fs shift into the surface
%     transmission mode (round edges won't change anything)
internalsettings.surf_tran_exp     =  2; % times wavelength @ fs transient zone width (cf. radar cross sections)
internalsettings.surf_tran_shift   =  0; % times wavelength @ fs shift out of the surface


% *******************************************************************************************************
% input parameter checks / prepare paramters

global globalsettings; % for sparse FIR performance characteristic

% number of input/output signals
ntx = size(settings, 1);
nrx = size(settings, 2);

% length of output signal => maximum length of all input signals
nout = max(cellfun(@length, input));
if nout == 0
   warn('No signals to transmit. Returning empty vectors (required length unknown).');
   output = cell(1, nrx);
   return;
end

% for all settings
for i = 1 : ntx
   for j = 1 : nrx
      % check contents of settings (substructs are needed => check everything)
      %     prepare required data
      expected.name = sprintf('settings{%.0f,%.0f}', i, j);
      expected.reqfields = {'probemode', 'type', 'c', 'fs',...
         'smallscale', {'on','det','seed','ensembles','maxrays','k','trms','eps_k','eps_trms','bw','fres','fs'},...
         'noise', {'type', 'n0', 'frxs', 'fs'}};
      %     either largescale+directivity or direct path setups (check for combination of both to be on the safe side)
      if isfield(settings{i,j}, 'direct') && ( isfield(settings{i,j}, 'largescale') || isfield(settings{i,j}, 'directivity') )
         err('Setup contains "direct" and "largescale/directivity" settings. This is not allowed.');
      elseif isfield(settings{i,j}, 'direct')
         expected.reqfields = [expected.reqfields, {'direct', {'gain', 'delay_s'}}];
      elseif ( isfield(settings{i,j}, 'largescale') || isfield(settings{i,j}, 'directivity') )
         expected.reqfields = [expected.reqfields, {...
            'largescale', {'type', 'dist', 'pl', 'c', 'f0', 'fs'},...
            'directivity', {'txant','txrot','txant','txrot','dir_tx'}}];
      end
      %     virtual connection (mirrored transmitter/receiver properties)
      expected.reqfields = [expected.reqfields, {'virtual', {'vtx','otx', 'vrx','orx', 'dim', 'refl'}}];
      %     reflective surfaces
      if isfield(settings{i,j}, 'surfaces')
         refl_names = fieldnames(settings{i,j}.surfaces);
         refl_cells = cell(0);
         for k = 1 : length(refl_names)
            refl_cells = [refl_cells, {refl_names{k}, {'mode', 'vtx', 'blur', 'n2', 'dim', 'pol_dim', 'aoi', 'poi2e', ...
               'largescale', {'type', 'dist', 'pl', 'c', 'f0', 'fs'},...
               'directivity', {'txant','txrot','txant','txrot','dir_tx'}}}];
         end
         expected.reqfields = [expected.reqfields, {'surfaces', refl_cells}];
         clear('refl_names', 'refl_cells');
      end
      %     check
      errortext = contentcheck(settings{i,j}, expected);
      %     output
      if ~isempty(errortext)
         err('Incomplete settings{%.0f,%.0f}\n%s', i, j, errortext);
      end
   end
end

% check if lines of vtx/otx and columns of vrx/orx are identical
vtx = cellfun(@(x) x.virtual.vtx, settings);
otx = cellfun(@(x) x.virtual.otx, settings);
vrx = cellfun(@(x) x.virtual.vrx, settings);
orx = cellfun(@(x) x.virtual.orx, settings);
if any( vtx(:,1) * size(vtx,2) ~= sum(vtx, 2) ) 
   err('At least one transmitter does not have identical .vtx settings for all receivers.');
end
if any( ( otx(:,1) * size(otx,2) ~= sum(otx, 2) ) & ~any(isnan(otx),2) ) 
   err('At least one transmitter does not have identical .otx settings for all receivers.');
end
if any( vrx(1,:) * size(vrx,1) ~= sum(vrx, 1) ) 
   err('At least one transmitter does not have identical .vrx settings for all receivers.');
end
if any( ( orx(1,:) * size(orx,1) ~= sum(orx, 1) ) & ~any(isnan(orx),1) ) 
   err('At least one transmitter does not have identical .orx settings for all receivers.');
end

% TX/RX should come first, then VTX/VRX; otherwise the function is inefficient
% (VTX: entire line is virtual; VRX: entire column is virtual)
isvtx = any( cellfun(@(x) x.virtual.vtx, settings), 2);
isvrx = any( cellfun(@(x) x.virtual.vrx, settings), 1);
if any(diff(isvtx) < 0)
   critwarn('Transmitters (TX) should come before virtual transmitters (VTX). This function is inefficient otherwise.');
end
if any(diff(isvrx) < 0)
   critwarn('Receivers (RX) should come before virtual receivers (VRX). This function is inefficient otherwise.');
end
ntx_nonv = sum(~isvtx);
nrx_nonv = sum(~isvrx);


% *******************************************************************************************************
% skip inactive transmitters by removing them and their VTX from the channel settings

% keep track of transmitter/receivers numbers, including VTX/VRX (temporary workaround until RX/TX are named)
txnum = 1 : ntx;
rxnum = 1 : nrx;

% for all real TX
for i = fliplr( find(isvtx==0)' ) % reversed order to make deleting entries possible
   % check if they do send a signal
   if ~isempty(input{i})
      continue 
   end 
   % find their mirror images
   otx = cellfun(@(x) x.virtual.otx, settings(:,1)); % lines are identical (check above)
   for k = fliplr( find(otx == i)' )
      settings(k,:) = [];
      txnum(k) = [];
   end
   settings(i,:) = [];
   txnum(i) = [];
end

% update number of transmitters
ntx = length(txnum);


% *******************************************************************************************************
% channel models: assemble channel impulse responses

for i = 1 : ntx        
   for j = 1 : nrx
      
      % vacuum these settings (speed up, remove "obsolete" settings, sanity check)
      settings{i,j} = vacuum(settings{i,j}, internalsettings, txnum(i), j);
      
      % large-scale (direct channel has to REPLACE largescale model + directivity; do not combine)
      if isfield(settings{i,j}, 'direct') % direct channel (results already provided)
         delay_ls = settings{i,j}.direct.delay_s;
         gain_ls  = settings{i,j}.direct.gain;
         if isfield(settings{i,j}.direct, 'gain_dir') % not mandatory; mostly created by vacuum() in this fcn
            gain_dir = settings{i,j}.direct.gain_dir;
         else
            gain_dir = 1;
         end
      else % normal model
         % largescale model
         [delay_ls, gain_ls] = channel_large(settings{i,j}.largescale);
         % directivity (antenna radiation patterns)
         gain_dir = channel_directivity(settings{i,j}.directivity);
      end 
      
      % surfaces (reflection)
      if isfield(settings{i,j}, 'surfaces')
         refl_names  = fieldnames(settings{i,j}.surfaces);
         % calculate delay and gain
         delay_rfl = zeros(length(refl_names), 1);
         gain_rfl  = zeros(length(refl_names), 1);
         for k = 1 : length(refl_names)
            [delay_rfl(k), gain_rfl(k)] = channel_large(settings{i,j}.surfaces.(refl_names{k}).largescale);
            gain_rfl(k) = gain_rfl(k) * channel_directivity(settings{i,j}.surfaces.(refl_names{k}).directivity);
            gain_rfl(k) = gain_rfl(k) * channel_surface(settings{i,j}.surfaces.(refl_names{k})) * ...
               settings{i,j}.surfaces.(refl_names{k}).gain;
            if settings{i,j}.surfaces.(refl_names{k}).mode == internalsettings.surf_tran
               criterr('Surface in transmission mode found in reflection block. Check nested function "vacuum".');
            end
         end
      end
      
      % small-scale
      if settings{i,j}.smallscale.on
         [delay_ss, gain_ss, statistics{i,j}] = channel_small(settings{i,j}.smallscale);
      else
         delay_ss = 0;
         gain_ss  = 1;
      end  
      
      % combine largescale and smallscale models
      gain{i,j} = gain_ls * gain_ss;
      delay_s{i,j} = delay_ls + delay_ss;
      %     directivity model
      switch lower(settings{i,j}.type)
         case 'room' % indoor, short-range: directivity applies only to LOS path
            gain{i,j}(1) = gain{i,j}(1) * gain_dir;
         case 'outdoor' % long-range: directivity applies to LOS and NLOS equally
            gain{i,j} = gain{i,j} * gain_dir;
         otherwise
            err('Unsupported settings.type "%s".', lower(settings.type));
      end
      
      % add surfaces (if present)
      if isfield(settings{i,j}, 'surfaces')
        [delay_s{i,j}, gain{i,j}] = merge_cirs(delay_s{i,j}, gain{i,j}, delay_rfl, gain_rfl);
      end
      
      % remove entries marked as obsolete
      delay_s{i,j}(isnan(gain{i,j})) = [];
      gain{i,j}(isnan(gain{i,j}))    = [];
      %     if not a single path survived: make sure the calculations below will still work 
      %      ... should be a very unlikely event anyway
      if isempty(delay_s{i,j})
         delay_s{i,j} = 1;
         gain{i,j}    = 0;
      end
      
      % user controlled gain factor
      if isfield(settings{i,j}, 'ucg')
         gain{i,j} = gain{i,j} * settings{i,j}.ucg;
      end
   end
end


% *******************************************************************************************************
% flatten CIRs by merging VTX and VRX

% which transmitters / receivers are virtual?
isvtx = cellfun(@(x) x.virtual.vtx, settings(:,1)); % identical for all RX (checks above)
isvrx = cellfun(@(x) x.virtual.vrx, settings(1,:)); % identical for all TX (checks above)
%     check if virtual transmitters and receivers are present at the same time (might indicate a setup problem)
if sum(isvtx) > 0 && sum(isvrx) > 1
   err('Both virtual transmitters and receivers present. This indicates a setup problem.'); % keep this an error for now
end

% get list of originating TX and RX
otx = cellfun(@(x) x.virtual.otx, settings(:,1)); % identical for all RX (checks above)
orx = cellfun(@(x) x.virtual.orx, settings(1,:)); % identical for all TX (checks above)

% flatten channels from real transmitters to real receivers
%     virtual transmitters
if sum(isvtx) > 0
   for i = find(isvtx==0)' % all real TX
      for j = 1 : nrx % all receivers
         for k = find(otx == i)' % all VTX created by this TX
            if sum(gain{k,j}) == 0; continue, end % skip dead connections
            [delay_s{i,j}, gain{i,j}] = merge_cirs(delay_s{i,j}, gain{i,j}, delay_s{k,j}, gain{k,j});
         end
      end
   end
end
%     virtual receivers
if sum(isvrx) > 0
   for j = find(isvrx==0) % all real RX
      for i = 1 : ntx % all transmitters
         for k = find(orx == j)% all VRX created by this RX
            if sum(gain{i,k}) == 0; continue, end % skip dead connections
            [delay_s{i,j}, gain{i,j}] = merge_cirs(delay_s{i,j}, gain{i,j}, delay_s{i,k}, gain{i,k});
         end
      end
   end
end
%     delete VTX and VRX from impulse responses and transmitter/receiver numbering
delay_s(isvtx,:) = [];
delay_s(:,isvrx) = [];
gain(isvtx,:) = [];
gain(:,isvrx) = [];
txnum(isvtx) = [];
rxnum(isvrx) = [];

% update number of remaining transmitters / receivers
ntx = length(txnum);
nrx = length(rxnum);


% *******************************************************************************************************
% record statistics

% initialize
statistics = cell(ntx_nonv, nrx_nonv);
statistics = cellfun(@(x) struct('pdp_los',NaN, 'pdp_sum',NaN, 'k',NaN, 'trms',NaN, 'tmin',NaN,...
   'tmax',NaN, 'cir_d',[], 'cir_g',[]), statistics, 'UniformOutput',false);

% record statistics (for active transmitters / receivers)
% % % close all; figure;
for i = 1 : ntx
   for j = 1 : nrx
% % %       subplot(ntx,nrx,(i-1)*nrx+j); stem(delay_s{i,j}, gain{i,j});
      statistics{txnum(i),rxnum(j)}.pdp_los = gain{i,j}(1).^2;
      statistics{txnum(i),rxnum(j)}.pdp_sum = sum( gain{i,j}.^2 );
      statistics{txnum(i),rxnum(j)}.k       = 10*log10( gain{i,j}(1)^2/sum(gain{i,j}(2:end).^2) );
      statistics{txnum(i),rxnum(j)}.trms    = sqrt(var((delay_s{i,j}-delay_ls)/settings{i,j}.smallscale.fs, gain{i,j}.^2));
      statistics{txnum(i),rxnum(j)}.tmin    = min(delay_s{i,j}) / settings{i,j}.smallscale.fs;
      statistics{txnum(i),rxnum(j)}.tmax    = max(delay_s{i,j}) / settings{i,j}.smallscale.fs;
      statistics{txnum(i),rxnum(j)}.cir_d   = delay_s{i,j} / settings{i,j}.smallscale.fs;
      statistics{txnum(i),rxnum(j)}.cir_g   = gain{i,j};
   end
end
% % % pause(1);
      
% if in probe mode: only return statistics; do not calculate filter
if any(any(cellfun(@(x) x.probemode, settings)))
   output = cell(1, nrx);
   return
end
      

% *******************************************************************************************************
% filter

% update number of transmittes and receivers
ntx = length(txnum);
nrx = length(rxnum);

% issue a warning if longest input signal is shorter than the transient region of the channel
len_inputs = cellfun(@length, input);
if max( cellfun(@(x) x.tmax/2, statistics) ) * settings{1,1}.smallscale.fs > min(len_inputs(len_inputs>0))
   critwarn('Transient region of channel is longer than the shortest input signal. The output signal will not reach stationarity!');
end

% output number of paths (only if at least one transmitter is active)
%     number of paths
if ntx > 1 || nrx > 1
   len = cellfun(@length, delay_s);
   msg('Considering the following amount of paths per channel transmitter -> receiver:');
   msg('        %s', sprintf('rx%02.0f ', [1:1:nrx]));
   for i = 1 : ntx
      if isempty(input{i}); continue; end % skip inactive transmitters
      msg('  tx%02.0f  %s', i, sprintf('%4.0f ', len(i,:)));
   end
else
   msg('Considering %.0f paths per channel transmitter -> receiver.', length(delay_s{1,1}));
end
%     warning for large amount of paths
if sum(sum(cellfun(@length, delay_s))) > internalsettings.maxpaths_warn
   warn('Large amount of paths to consider; this may take a while.');
end

% load performance characteristic to select fastest filtering method
%  ... make sure "stubborn matfile loading" (ignore missing files) is enabled
mat_stubborn_temp = globalsettings.core.mat_stubborn;
globalsettings.core.mat_stubborn = true;
pchar = loadmat(globalsettings.core.sfir_char);
globalsettings.core.mat_stubborn = mat_stubborn_temp;
pchar.exists = ~isempty(fieldnames(pchar));
if ~pchar.exists % no performance char file exists for this host
   warn('No sparse FIR performance found for this host. Switching to defaults.');
   pchar.xing_pt = internalsettings.default_xing_pt;
end

% apply delay and gain
for j = 1 : nrx
   output{j} = zeros(nout, 1);
   for i = 1 : ntx
      % input signal length
      len = length(input{txnum(i)});
      if len == 0; continue; end % skip inactive transmitters
      
      % get performance crossing between filter and loop implementation (if available)
      if pchar.exists
         ind_i = interp1(pchar.len_i_log10, pchar.ind_i, log10(len), 'nearest', 'extrap');
         ind_f = interp1(pchar.len_f_log10, pchar.ind_f, log10(delay_s{i,j}(end)), 'nearest', 'extrap');
         pchar.xing_pt = pchar.sxing(ind_i, ind_f);
      end
      
      % filtering
      if (delay_s{i,j}(end)-length(delay_s{i,j}))/delay_s{i,j}(end) > pchar.xing_pt % loop faster (low density)
         for k = 1 : length(delay_s{i,j}) % for all paths (assumed to be sparse => didn't use filter)
            output{j}(1+delay_s{i,j}(k):len) =...
               output{j}(1+delay_s{i,j}(k):len) + input{txnum(i)}(1:end-delay_s{i,j}(k)) * gain{i,j}(k);
         end
      else
         gain_full = zeros(1+delay_s{i,j}(end), 1); % first term in numerator of FIR corresponds to delay_s=0
         gain_full(1+delay_s{i,j}) = gain{i,j};
         output{j}(1:len) = output{j}(1:len) + filter(gain_full, 1, input{txnum(i)});
      end
   end
end


% *******************************************************************************************************
% noise
%  ... using cell2mat / mat2cell and the ability of channel_noise to process multiple receivers at
%      once is slower for small amounts of receivers => treat all receivers separately

for j = 1 : nrx
   if ~strcmpi(settings{j}.noise.type, 'off')
      output{j} = channel_noise(output{j}, settings{1,j}.noise); % identical for all transmitters
   end
end

end



% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% vacuum settings (remove contradictory or unnecessary parts plus do some sanity checks)
% ... settings are not returned by channel_main => no external effects of these changes
% NOTE THAT the treatment of virtual transmitters in here is an approximation and works only for
%           parallel surfaces. This part should be replaced by automatically generated VTX using
%           raytracing.
function settings = vacuum(settings, internalsettings, tx, rx)

% extract some information (speed up checks below)
%     booleans
internal.surfaces = isfield(settings, 'surfaces');
internal.virtual  = settings.virtual.vtx;
%     surface names (if present)
if internal.surfaces
   internal.surf_names = fieldnames(settings.surfaces);
end

% reflective surfaces without a propagation model (this should not happen)
if internal.surfaces && ~(isfield(settings, 'largescale') && isfield(settings, 'directivity'))
   critwarn('The channel has reflective surfaces, but no largescale/directivity model; removing the surfaces.');
   settings = rmfield(settings, 'surfaces');
end

% calculate rolloff factor for surfaces (additional gain 0..1 that ensures a smooth transistion between
% modes "transmit"/"reflect" and "off" => simple diffraction/refraction model)
if internal.surfaces
   for i = 1 : length(internal.surf_names)      
      % POI to closest edge in wavelengths at carrier level
      poi2e_lambda = settings.surfaces.(internal.surf_names{i}).poi2e / settings.c * settings.fs;
      % smooth rolloff if POI is close to the edge ("diffraction")
      if settings.surfaces.(internal.surf_names{i}).mode == internalsettings.surf_tran
         % transmission: diffraction will cause wave prop. behind the surface (=> -poi2e_lambda)
         settings.surfaces.(internal.surf_names{i}).gain = ...
            exp(- max(0, (-poi2e_lambda + internalsettings.surf_tran_shift * settings.surfaces.(internal.surf_names{i}).blur) /...
            (internalsettings.surf_tran_exp * settings.surfaces.(internal.surf_names{i}).blur) ));
      else
         % reflection: diffraction will cause wave prop. outside the surface's immediate refl. area
         settings.surfaces.(internal.surf_names{i}).gain = ...
            exp(- max(0, (poi2e_lambda + internalsettings.surf_refl_shift * settings.surfaces.(internal.surf_names{i}).blur) /...
            (internalsettings.surf_refl_exp * settings.surfaces.(internal.surf_names{i}).blur) ));
      end
   end
end

% virtual connection
%    1a) switch off all surfaces that have a "covered by VTX" tag and are in reflect mode
%    1b) create a warning if there are surfaces in the spatial dimension of the VTX that are not tagged
%        (avoid reflection/transmission by surfaces that created the VTX: "mirror of mirror = original")
%    2a) if there is no linked surface: connection always active
%    2b) if there is a linked surface: kill connection if does not pass through (i.e., surf. in transmit mode)
if internal.virtual
   % there is a linked surface
   if isfield(settings.virtual, 'surface')
      % surface not found
      if ~internal.surfaces || ~isfield(settings.surfaces, settings.virtual.surface)
         err('Linked surface "%s" for virtual transmitter %i not found.', settings.virtual.surface, tx);
      end
      % "inside the mirror" (this can happen for doors, but can also be a setup problem)
      if settings.surfaces.(settings.virtual.surface).mode == internalsettings.surf_refl
         warn('Virtual transmitter %i and receiver %i are on the same side of the mirror surface "%s". Removing the VTX.',...
            tx, rx, settings.virtual.surface);
         % remove the connection
         kill_connection();
         return
      end
      % should the linked surface be inverted?
      if isfield(settings.virtual, 'invsurf') && settings.virtual.invsurf
         settings.surfaces.(settings.virtual.surface).gain = sqrt(1 - settings.surfaces.(settings.virtual.surface).gain.^2);
      end
      % if the link VTX -> RX passes through the linked surface's plane (surface is in transmit mode)
      %     yes: copy surface gain factor to the link gain and switch it off
      if settings.surfaces.(settings.virtual.surface).mode == internalsettings.surf_tran
         settings.ucg = settings.ucg * settings.surfaces.(settings.virtual.surface).gain;
         settings.surfaces.(settings.virtual.surface).mode = internalsettings.surf_off;
      %     no: kill the connection (including surfaces)
      else
         kill_connection();
         return
      end
   end
   
   % switch off all surfaces that have been tagged "covered by VTX" and check for "dangerous" setups
   for i = 1 : length(internal.surf_names)
      if settings.surfaces.(internal.surf_names{i}).vtx
         settings.surfaces.(internal.surf_names{i}).mode = internalsettings.surf_off;
      % reflecting surface left in the dimension of the VTX: most likely a setup problem
      elseif settings.surfaces.(internal.surf_names{i}).dim == settings.virtual.dim &&...
            settings.surfaces.(internal.surf_names{i}).mode == internalsettings.surf_refl
         critwarn(['Surface "%" in dim of VTX %i is not marked as "covered by VTX" and set to reflection.', ...
            'This might be a setup problem.'], internal.surf_names{i}, tx);
      end
   end
end

% switch off all surfaces with a gain below the threshold and get modes of surfaces
if internal.surfaces
   for i = 1 : length(internal.surf_names)      
      if settings.surfaces.(internal.surf_names{i}).mode == internalsettings.surf_refl && ...
            settings.surfaces.(internal.surf_names{i}).gain < internalsettings.surf_mingain
         settings.surfaces.(internal.surf_names{i}).gain = 0;
         settings.surfaces.(internal.surf_names{i}).mode = internalsettings.surf_off;
      elseif settings.surfaces.(internal.surf_names{i}).mode == internalsettings.surf_tran && ...
            settings.surfaces.(internal.surf_names{i}).gain > (1 - internalsettings.surf_mingain)
         settings.surfaces.(internal.surf_names{i}).gain = 1;
         settings.surfaces.(internal.surf_names{i}).mode = internalsettings.surf_off;
      end
      internal.surf_mode(i) = settings.surfaces.(internal.surf_names{i}).mode;
   end
end

% if there are surfaces in "transmit mode" => direct path blocked by at least one surface
if internal.surfaces && sum(internal.surf_mode == internalsettings.surf_tran) > 0
   % replace largescale model by sum of all transmitting surfaces
   %     get LOS channel (gain and delay)
   [delay_s, gain] = channel_large(settings.largescale);
   %     VTX: surface is passed several times (2,2,4,4,... for 1,2,..-pt refl)
   if internal.virtual
      surf_exp = ceil(settings.virtual.refl / 2) * 2;
   else
      surf_exp = 1;
   end
   %     get gains of all transmitting surfaces and switch them off
   for i = find(internal.surf_mode == internalsettings.surf_tran)      
      gain = gain * ( settings.surfaces.(internal.surf_names{i}).gain + ...
         (1 - settings.surfaces.(internal.surf_names{i}).gain) * ...
         channel_surface(settings.surfaces.(internal.surf_names{i})) ^ surf_exp);
      settings.surfaces.(internal.surf_names{i}).mode = internalsettings.surf_off;
   end  
   %     set this gain/delay in "direct mode"
   settings.direct.gain     = gain;
   settings.direct.gain_dir = channel_directivity(settings.directivity);
   settings.direct.delay_s  = delay_s;
   %     and remove largescale and directivity model (keep smallscale!)
   settings = rmfield(settings, {'largescale', 'directivity'});
end

% remove inactive surfaces (will not contribute anyway)
if internal.surfaces
   for i = 1 : length(internal.surf_names)
      if settings.surfaces.(internal.surf_names{i}).mode == internalsettings.surf_off;
         settings.surfaces     = rmfield(settings.surfaces, internal.surf_names{i});
         internal.surf_mode(i) = internalsettings.surf_off;
      end
   end
   internal.surf_names(internal.surf_mode == internalsettings.surf_off) = [];
   % if there are no surfaces left: remove the struct altogether
   if isempty(fieldnames(settings.surfaces))
      settings = rmfield(settings, 'surfaces');
   end
end


   % ****************************************************************************************************
   % remove a connection completely
   function kill_connection
      settings.direct.gain     = 0;
      settings.direct.gain_dir = 0;
      settings.direct.delay_s  = 1;
      settings.smallscale.on   = false;
      settings = rmfield(settings, {'largescale', 'directivity'});
      if isfield(settings, 'surfaces')
         settings = rmfield(settings, 'surfaces');
         internal.surf_names = {};
         internal.surfaces   = false;
      end
   end
end
   


% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% merge two channel impulse responses
function [delay, gain] = merge_cirs(delay1, gain1, delay2, gain2)

% insert the shorter CIR into the longer
% (use insertion sort, since delay bins might already be present)
if length(delay1) > length(delay2)
   [delay, gain] = insertion_sort(delay1, gain1, delay2, gain2);
else
   [delay, gain] = insertion_sort(delay2, gain2, delay1, gain1);
end

  % merge gi(di) into d(g) using insertion sort
  function [d, g] = insertion_sort(d, g, di, gi)
    % insert/merge
    for k = 1 : length(di)
      ind_k = find(d == di(k));
      switch length(ind_k)
        case 0 % not found => just add and sort later in one operation
          d = [d; di(k)];
          g = [g; gi(k)];
        case 1 % found => add gains directly (allows for constructive and destructive interference)
           g(ind_k) = g(ind_k)^2 + gi(k)^2;
%           g(ind_k) = sqrt(g(ind_k)^2 + gi(k)^2) * sign(g(ind_k)) * sign(gi(k)); % add power
        otherwise % found several entries => add to the first one
           g(ind_k(1)) = g(ind_k(1))^2 + gi(k)^2;
%           g(ind_k(1)) = sqrt(g(ind_k(1))^2 + gi(k)^2) * sign(ind_k(1)) * sign(gi(k)); % add power
      end
    end
    % resort
    [d, sort_ind] = sort(d);
    g = g(sort_ind);
  end
end


            
            
% *******************************************************************************************************
% *******************************************************************************************************
% OLD AND RUSTY

%       % add surfaces (if present)
%       if isfield(settings{i,j}, 'surfaces')
%          % use insertion sort (delay bins might already be present)
%          for k = 1 : length(delay_rfl)
%             ind_k = find(delay_s{i,j} == delay_rfl(k));
%             switch length(ind_k)
%                case 0 % not found => just add and sort later
%                   delay_s{i,j} = [delay_s{i,j}; delay_rfl(k)];
%                   gain{i,j}    = [gain{i,j}; gain_rfl(k)];
%                case 1 % found => add power
%                   gain{i,j}(ind_k) = sqrt(gain{i,j}(ind_k)^2 + gain_rfl(k)^2);
%                otherwise % shouldn't happen, but might pose a performance problem (no error though!)
%                   critwarn('Delay vector contains identical entries (TX%i, RX%i).', i,j);
%                   gain{i,j}(ind_k(1)) = sqrt(gain{i,j}(ind_k(1))^2 + gain_rfl(k)^2);
%             end
%          end
%          %     resort
%          [delay_s{i,j}, sort_ind] = sort(delay_s{i,j});
%          gain{i,j} = gain{i,j}(sort_ind);
%       end

%       % add surfaces (if present)
%       if isfield(settings{i,j}, 'surfaces')        
%          delay_s{i,j} = [delay_s{i,j}; delay_rfl];
%          gain{i,j}    = [gain{i,j}; gain_rfl];
%          % sort delays
%          [delay_s{i,j}, sort_ind] = sort(delay_s{i,j});
%          gain{i,j} = gain{i,j}(sort_ind);
%          % remove entries marked as obsolete
%          delay_s{i,j}(isnan(gain{i,j})) = [];
%          gain{i,j}(isnan(gain{i,j}))    = [];
%          % make sure there are no identical delays
%          if any(diff(delay_s{i,j})==0)
%             for k = 1 : length(delay_s{i,j}) - 1 % [!!!] FIXME
%                if delay_s{i,j}(k) == delay_s{i,j}(k+1)
%                   gain{i,j}(k+1) = sqrt(gain{i,j}(k+1)^2 + gain{i,j}(k)^2);
%                   gain{i,j}(k)   = NaN; % mark as obsolete
%                end
%             end
%          end
%       end

% % %       settings.surfaces.(internal.surf_names{i}).gain = ...
% % %          1/2 + 1/2 * cos(pi * (min(1, max(0,
% (poi2e_lambda+internalsettings.surf_refl_shift)/internalsettings.surf_refl_exp))) );

%    % replace largescale model by linked surface
%    %     get channel (gain and delay) including directivity
%    [delay_s, gain] = channel_large(settings.surfaces.(settings.virtual.surface).largescale);
%    gain = gain * channel_directivity(settings.surfaces.(settings.virtual.surface).directivity);
%    gain = gain * channel_surface(settings.surfaces.(settings.virtual.surface));
%    %     set this gain/delay in "direct mode"
%    settings.direct.gain    = gain;
%    settings.direct.delay_s = delay_s;
%    %     remove largescale and directivity model (keep smallscale!)
%    settings = rmfield(settings, {'largescale', 'directivity'});
%    %     and also the linked surface (now in direct channel model)
%    settings.surfaces = rmfield(settings.surfaces, settings.virtual.surface);

%             % find identical entries [start1, stop1, start2, stop2, ...]
%             ident_bnds = findzeros( [1; diff(delay_s{i,j}); 1] - 0.1 ); % -0.1 => bias towards identical
%             % sum up / construct an index array for swift removal of all identical entries
%             ident_ind = zeros(sum(ident_bnds(2:2:end)-ident_bnds(1:2:end-1)+1), 1);
%             ident_ptr = 1;
%             for k = 1 : length(ident_bnds)/2
%                gain{i,j}(ident_bnds(2*k-1)-1) = sum(gain{i,j}(ident_bnds(2*k-1)-1:ident_bnds(2*k)));
%                ident_ind(ident_ptr:ident_ptr+ident_bnds(2*k)-ident_bnds(2*k-1)) = ident_bnds(2*k-1):ident_bnds(2*k);
%                ident_ptr = ident_ptr + ident_bnds(2*k) - ident_bnds(2*k-1) + 1;
%             end
%             % clean up
%             delay_s{i,j}(ident_ind) = [];
%             gain{i,j}(ident_ind) = [];
%             clear('ident_bnds', 'ident_ind', 'ident_ptr');

