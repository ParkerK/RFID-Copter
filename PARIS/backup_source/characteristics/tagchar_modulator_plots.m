% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate plots for: tag - modulation (reflection coefficient)
%
% WARNING: Settings are not checked for sanity!
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
clear; close all; clc; pause(0.01);
version = 'beta 3.0';

% paths, etc
%     path to globalinit.m
path = fullfile(fileparts(mfilename('fullpath')), '..');
%     Matlab does not resolve symbolic links; if we're running on Unix, let the system resolve them
if isunix
   [dummy, path] = system(sprintf('readlink -f "%s"', path));
end
%     add this path and call globalinit script
restoredefaultpath; addpath(path); clear('dummy', 'path');
globalinit('silent');

% "switch to here"
cd(fileparts(mfilename('fullpath')));

% output own name/version
disp(sprintf('*******************************************************************************************************'));
disp(sprintf('* This is %s, version %s', mfilename, version));
disp(sprintf('*******************************************************************************************************'));


% *******************************************************************************************************
% setup

% assembly and detuning state mapping (filename)
settings.charfiles.ads_mapping = 'tagchar_modulator_adsm';

% characteristic ID (prefix)
settings.charid     = '_c-p2a'; % start with '_' if not empty

% assembly and detuning states
cat = 450e-15; %[450, 1250]; % F
rat =       0; %[-67, 0, 200]; % percent shift
enr = 0;
fsr = 0;
% enr =     0.5; %[0, 0.5]; % boost of resonance (slightly weakened res. up to near-metal)
% fsr =   100e6; %[0, 100e6]; % frequency shift of resonance

% map assembly and detuning to states
%     load mapping
ads_mapping = load(settings.charfiles.ads_mapping);
%     map
[settings.state_a, settings.state_d] = adstate('settings->state', ads_mapping, cat, rat, enr, fsr);
% settings.state_a = find( round(ads_mapping.assembly.v_cat*1e15) == cat & ads_mapping.assembly.v_rats == rat );
% settings.state_d = find( abs(ads_mapping.detuning.v_res_en - enr) < 1e-10 & ads_mapping.detuning.v_fshift == fsr );
%     check
if isempty(settings.state_a)
   error('Assembly impedance not part of the state vectors.');
end
if isempty(settings.state_d)
   error('Detuning state not part of the state vectors.');
end

% load workspace of tagchar_modulator
%     complete filename
settings.suffix = sprintf('%s-sim_a%03i_d%03i', settings.charid, settings.state_a, settings.state_d);
%     backup own version number ("version" also in mat-file)
ownversion = version;
%     initialize f as vector ... otherwise f() will lead to an error
f = [];
load(sprintf('tagchar_mod%s_workspace', settings.suffix));
%     redefine rho and rho_pav
rho_pin = rho;
rho     = rho_pav;
clear rho_pav;
%     version stuff
version = ownversion; clear ownversion;

% print what is done
disp(sprintf('\n\n*******************************************************************************************************'));
disp(sprintf('* WORKING ON: %s.mat (created by %s)', matfilename, createdby));

% plot position of nwarperr largest warping errors
settings.nwarperr = 100;

% minimum/maximum for some plots
settings.f_min = 0.50; % GHz
settings.f_max = 1.30; % GHz
settings.pav_min = 1e-6; % W
settings.pav_max =    3; % W

% erosion of border region (e.g. for drho/df) => eroded reflection coefficient
settings.erode = 3; % values from border

% clim for some plots
%     have to be set
settings.clim.delta_rho = [0, 2/3]; % |delta_rho| vs f and pav (xy-plane)
settings.clim.rho_c_mag = [0,   1]; % rho_c vs f and pav (linear model)
settings.clim.rho_d_mag = [0,   1]; % rho_d vs f and pav (linear model)
settings.clim.rho_c_arg = [   0, 180]; % deg rho_c vs f and pav (linear model)
settings.clim.rho_d_arg = [-180, 180]; % deg rho_d vs f and pav (linear model)
%     can be empty ([] for auto)
settings.clim.reserr         = []; % residual IIR fitting error [-50,-15]
settings.clim.mfcw_derr1     = [ -1,  1]; % m 2FCW distance error (full span plot)
settings.clim.mfcw_derr2     = []; % cm 2FCW distance error (zoomed plot)
settings.clim.mfcw_derr_hist = [0, 30]; % percent 2FCW distance error histogram
settings.clim.t_fitiir       = []; % s IIR fitting time
settings.clim.darg_rho_d1    = []; % deg zoomed phase difference: rho_d - rho_d(settings.f1)
settings.clim.darg_rho_d2    = []; % deg zoomed phase difference: rho_d - rho_d(settings.f2)
settings.clim.darg_rho_d3    = []; % deg zoomed phase difference: rho_d - rho_d(settings.f3)
settings.clim.darg_rho_d4    = []; % deg zoomed phase difference: rho_d - rho_d(settings.f4)

% tick for some plots (where this is necessary)
settings.pav_tick  = [1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0]; % W
settings.fop_tick  = settings.f_min : 0.1 : settings.f_max; % GHz

% pav for example filter characteristic plots
settings.pav1 = 2e-4; % W
settings.pav2 = 1e-2; % W

% assumed carrier center frequencies (for some plots)
settings.f1 = 865.7e6; % Hz (European frequency band)
settings.f2 = 867.5e6; % Hz (European frequency band)
settings.f3 = 902.0e6; % Hz (US frequency band)
settings.f4 = 928.0e6; % Hz (US frequency band)

% plot secondary carrier @ settings.f1 to .f4 +/- settings.fd (for some linear model plots)
settings.fd = 50e6; % Hz

% MFCW distance error estimates
settings.mfcw.df12        = [1e6, 10e6]; % Hz frequency differences between 2FCW carrier pairs (plot: first two)
settings.mfcw.hist_n      =          50; % number of points for distance error histogram
settings.mfcw.hist_minmax = [-1.5, 1.5]; % m bounds for distance error histogram

% output folder for plots
settings.plotfolder = 'plots'; % subfolder for plots

% figure options
settings.figopt.papersize       = [21, 12.5]; % cm size [width, height] 
settings.figopt.visible         = 'off'; % figure visibility {'on', 'off'} 
settings.figopt.type            = 'png'; % type of figures to create ('png' or 'eps'; leave empty to create both; see SAVEFIGURE)
settings.figopt.nonconfidential = false; % remove confidential information (e.g. for presentations)?
settings.figopt.resolution      =   200; % dpi image resolution
%     settings for thesis
settings.figopt.papersize       = [14.5, 10]; % cm size [width, height] 
settings.figopt.visible         = 'off'; % figure visibility {'on', 'off'} 
settings.figopt.type            = 'eps'; % type of figures to create ('png' or 'eps'; leave empty to create both; see SAVEFIGURE)
settings.figopt.nonconfidential = false; % remove confidential information (e.g. for presentations)?
settings.figopt.resolution      =   200; % dpi image resolution
settings.figopt.addspace_ylabel =  0.00;
settings.figopt.addspace_colorb =  0.03;




% *******************************************************************************************************
% prepare data: indices

ind.fopmin = interp1(fch, [1:1:settings.nfch], settings.fopmin, 'nearest', 'extrap');
ind.fopmax = interp1(fch, [1:1:settings.nfch], settings.fopmax, 'nearest', 'extrap');
ind.pav1   = interp1(pav, [1:1:settings.np],   settings.pav1,   'nearest', 'extrap');
ind.pav2   = interp1(pav, [1:1:settings.np],   settings.pav2,   'nearest', 'extrap');
ind.f1     = interp1(f,   [1:1:settings.nf],   settings.f1,     'nearest', 'extrap');
ind.fch1   = interp1(fch, [1:1:settings.nfch], settings.f1,     'nearest', 'extrap');
ind.fch2   = interp1(fch, [1:1:settings.nfch], settings.f2,     'nearest', 'extrap');
ind.fch3   = interp1(fch, [1:1:settings.nfch], settings.f3,     'nearest', 'extrap');
ind.fch4   = interp1(fch, [1:1:settings.nfch], settings.f4,     'nearest', 'extrap');


% *******************************************************************************************************
% prepare data: eroded reflection coefficient
disp(sprintf('\nEroding reflection coefficient borders...'));

% erode reflection coefficient borders to remove high gradients at borders
%     prepare matrices and copy data (add NaNs to borders)
rhoe = nan(settings.nfch+2*settings.erode, settings.nm, settings.np+2*settings.erode);
rhoe(1+settings.erode:settings.erode+settings.nfch, :, 1+settings.erode:settings.erode+settings.np) = rho;
eroded_f = rhoe;
eroded_p = rhoe;
%     erode ("filter" with window of size 2*settings.erode+1, set center to NaN if any NaN present)
for i = 1 : settings.nm % work in slices for zmod (does not have to be eroded)
   % erode along fch
   for j = 1+settings.erode : settings.np+settings.erode
      for k = 1+settings.erode : settings.nfch+settings.erode;
         if any(isnan( rhoe(k-settings.erode:k+settings.erode, i, j) ))
            eroded_f(k, i, j) = NaN;
         end
      end
   end
   % erode along p
   for j = 1+settings.erode : settings.nfch+settings.erode
      for k = 1+settings.erode : settings.np+settings.erode;
         if any(isnan( rhoe(j, i, k-settings.erode:k+settings.erode) ))
            eroded_p(j, i, k) = NaN;
         end
      end
   end
end
%     combine (_f and _p are either equal or NaN)
rhoe = (eroded_f + eroded_p) / 2;
%     extract data, compensate for shift
rhoe = rhoe(1+settings.erode:settings.erode+settings.nfch, :, 1+settings.erode:settings.erode+settings.np);

% cleanup
clear('eroded_f', 'eroded_p');


% *******************************************************************************************************
% prepare data: linear model (rho_c + rho_{delta})
disp(sprintf('Preparing data for linear model...'));

% reflection coefficient
rho_c = ( squeeze(rho(:, 1, :)) + squeeze(rho(:, end, :)) ) / 2;
rho_d =   squeeze(rho(:, 1, :)) - rho_c;

% eroded reflection coefficient
rhoe_c = ( squeeze(rhoe(:, 1, :)) + squeeze(rhoe(:, end, :)) ) / 2;
rhoe_d =   squeeze(rhoe(:, 1, :)) - rhoe_c;


% *******************************************************************************************************
% prepare data: reflection coefficient magnitude / phase
disp(sprintf('Preparing reflection coefficient magnitude / phase...'));

% magnitue
mag.rho     = abs(rho);
mag.rho_ext = abs(rhoext_pav);
mag.rhoe    = abs(rhoe);
mag.rho_c   = abs(rho_c);
mag.rho_d   = abs(rho_d);
mag.rhoe_c  = abs(rhoe_c);
mag.rhoe_d  = abs(rhoe_d);

% phase
arg.rho     = angle(rho);
arg.rho_ext = angle(rhoext_pav);
arg.rhoe    = angle(rhoe);
arg.rho_c   = angle(rho_c);
arg.rho_d   = angle(rho_d);
arg.rhoe_c  = angle(rhoe_c);
arg.rhoe_d  = angle(rhoe_d);
%     unwrap
for i = [2,3,1]
   arg.rho    (end:-1:1,:,:) = unwrap(arg.rho    (end:-1:1,:,:), [], i);
   arg.rho_ext(end:-1:1,:,:) = unwrap(arg.rho_ext(end:-1:1,:,:), [], i);
   arg.rhoe   (end:-1:1,:,:) = unwrap(arg.rhoe   (end:-1:1,:,:), [], i);
   arg.rho_c  (end:-1:1,:,:) = unwrap(arg.rho_c  (end:-1:1,:,:), [], i);
   arg.rho_d  (end:-1:1,:,:) = unwrap(arg.rho_d  (end:-1:1,:,:), [], i);
   arg.rhoe_c (end:-1:1,:,:) = unwrap(arg.rhoe_c (end:-1:1,:,:), [], i);
   arg.rhoe_d (end:-1:1,:,:) = unwrap(arg.rhoe_d (end:-1:1,:,:), [], i);
end
%     rad -> deg
arg.rho     = arg.rho     * (180/pi);
arg.rho_ext = arg.rho_ext * (180/pi);
arg.rhoe    = arg.rhoe    * (180/pi);
arg.rho_c   = arg.rho_c   * (180/pi);
arg.rho_d   = arg.rho_d   * (180/pi);
arg.rhoe_c  = arg.rhoe_c  * (180/pi);
arg.rhoe_d  = arg.rhoe_d  * (180/pi);


% *******************************************************************************************************
% prepare data: reflection coefficient gradient drho / df (use eroded reflection coefficient)
disp(sprintf('Calculating reflection coefficient magnitute/phase gradients...'));

% create fch vector for drho (center values)
fch_drho = (fch(1:end-1) + fch(2:end)) / 2;

% df
dfch = diff(fch);

% reflection coefficient
%     prepare matrices
mag.drho = nan(settings.nfch-1, settings.nm, settings.np);
arg.drho = nan(settings.nfch-1, settings.nm, settings.np);
%     calculate d/df
for i = 1 : settings.nm
   for j = 1 : settings.np
      mag.drho(:, i, j) = diff(mag.rho(:, i, j)) ./ dfch;
      arg.drho(:, i, j) = diff(arg.rho(:, i, j)) ./ dfch;
   end
end

% eroded reflection coefficient
%     prepare matrices
mag.drhoe = nan(settings.nfch-1, settings.nm, settings.np);
arg.drhoe = nan(settings.nfch-1, settings.nm, settings.np);
%     calculate d/df
for i = 1 : settings.nm
   for j = 1 : settings.np
      mag.drhoe(:, i, j) = diff(mag.rhoe(:, i, j)) ./ dfch;
      arg.drhoe(:, i, j) = diff(arg.rhoe(:, i, j)) ./ dfch;
   end
end

% linear model
%     prepare matrices
mag.drho_c = nan(settings.nfch-1, settings.np);
arg.drho_c = nan(settings.nfch-1, settings.np);
mag.drho_d = nan(settings.nfch-1, settings.np);
arg.drho_d = nan(settings.nfch-1, settings.np);
%     calculate d/df
for j = 1 : settings.np
   mag.drho_c(:, j) = diff(mag.rho_c(:, j)) ./ dfch;
   arg.drho_c(:, j) = diff(arg.rho_c(:, j)) ./ dfch;
   mag.drho_d(:, j) = diff(mag.rho_d(:, j)) ./ dfch;
   arg.drho_d(:, j) = diff(arg.rho_d(:, j)) ./ dfch;
end

% linear model, eroded reflection coefficient
%     prepare matrices
mag.drhoe_c = nan(settings.nfch-1, settings.np);
arg.drhoe_c = nan(settings.nfch-1, settings.np);
mag.drhoe_d = nan(settings.nfch-1, settings.np);
arg.drhoe_d = nan(settings.nfch-1, settings.np);
%     calculate d/df
for j = 1 : settings.np
   mag.drhoe_c(:, j) = diff(mag.rhoe_c(:, j)) ./ dfch;
   arg.drhoe_c(:, j) = diff(arg.rhoe_c(:, j)) ./ dfch;
   mag.drhoe_d(:, j) = diff(mag.rhoe_d(:, j)) ./ dfch;
   arg.drhoe_d(:, j) = diff(arg.rhoe_d(:, j)) ./ dfch;
end


% *******************************************************************************************************
% difference of phase to some center value (for linear model)
disp(sprintf('Phase/Magnitude difference to carriers (linear model) ...'));

% carrier at settings.f1
linmodel.darg.rho_c1 = arg.rho_c - repmat(arg.rho_c(ind.fch1, :), size(fch), 1);
linmodel.darg.rho_d1 = arg.rho_d - repmat(arg.rho_d(ind.fch1, :), size(fch), 1);
% carrier at settings.f2
linmodel.darg.rho_c2 = arg.rho_c - repmat(arg.rho_c(ind.fch2, :), size(fch), 1);
linmodel.darg.rho_d2 = arg.rho_d - repmat(arg.rho_d(ind.fch2, :), size(fch), 1);
% carrier at settings.f3
linmodel.darg.rho_c3 = arg.rho_c - repmat(arg.rho_c(ind.fch3, :), size(fch), 1);
linmodel.darg.rho_d3 = arg.rho_d - repmat(arg.rho_d(ind.fch3, :), size(fch), 1);
% carrier at settings.f4
linmodel.darg.rho_c4 = arg.rho_c - repmat(arg.rho_c(ind.fch4, :), size(fch), 1);
linmodel.darg.rho_d4 = arg.rho_d - repmat(arg.rho_d(ind.fch4, :), size(fch), 1);


% *******************************************************************************************************
% prepare data: minimum available power levels
disp(sprintf('Calculating minimum available power necessary for tag operation...'));

pavopmin.all      = zeros(settings.nfch, settings.nm);
pavopmin.square50 = zeros(settings.nfch,           1);
for i = 1 : settings.nfch
   
   % Pav,min vs. frequency and state of modulation
   for j = 1 : settings.nm
      pavopmin.all(i,j) = interp1(pic, squeeze(m_pav(i, j, :)), picopmin(i), 'linear', 'extrap'); 
   end
   
   % Pav,min vs. frequency for perfect square modulation w. 50% duty cycle
   for j = 1 : settings.np
      %     pav -> pic for [mod, unmod]
      pic10(j, 1) = interp1(squeeze(m_pav(i,   1, :)), pic, pav(j), 'linear', 'extrap');
      pic10(j, 2) = interp1(squeeze(m_pav(i, end, :)), pic, pav(j), 'linear', 'extrap');
   end
   %     average pic -> pav
   pavopmin.square50(i) = interp1(mean(pic10,2), pav, picopmin(i), 'spline', 'extrap');
end


% *******************************************************************************************************
% prepare data: background for plots (higher values for more support, NaN: no support)
disp(sprintf('Preparing background for plots...'));

% for (fch, pav) plots
background.f_pav = ones(settings.nfch, settings.np) * settings.nm;
for i = 1 : settings.nm
   background.f_pav = background.f_pav - isnan(squeeze(mag.rho(:, i, :))); % 0 ... settings.nm
end
%     create a version with NaN outside support
background.f_pav_nan = background.f_pav;
background.f_pav_nan(find(background.f_pav==0)) = NaN; %#ok<FNDSB>
%     normalize -1...0
background.f_pav_nan = background.f_pav_nan / max(max(background.f_pav)) - 1;
background.f_pav     = background.f_pav     / max(max(background.f_pav)) - 1;


% *******************************************************************************************************
% prepare data: warping error
disp(sprintf('Processing warping error data...'));

% sort error data by average (descending)
[err_sorted, sort_ind] = sort(abs(err(:)), 'descend');
%   align indices
err_ind = repmat(err_ind, 2, 1); % not 100 largest "real/imag part errors", but largest "overall errors"
ind.err = err_ind(sort_ind, :);

% z-data for error positions (in mag/phase plots)
mag.rhoerr = nan(settings.nwarperr, 2); % ["all", >maxerr]
arg.rhoerr = nan(settings.nwarperr, 2); % ["all", >maxerr]
for i = 1 : settings.nwarperr
    % "all"
    mag.rhoerr(i, 1) = mag.rho(ind.err(i,1), ind.err(i,2), ind.err(i,3));
    arg.rhoerr(i, 1) = arg.rho(ind.err(i,1), ind.err(i,2), ind.err(i,3));
    % > maxerr
    if err_sorted(i) > settings.maxerr
       mag.rhoerr(i, 2) = mag.rhoerr(i, 1);
       arg.rhoerr(i, 2) = arg.rhoerr(i, 1);
    end
end


% *******************************************************************************************************
% MFCW distance errors based on linear model
% d_err = (phi_d0-phi_d1) * c/(2*w_d);
disp(sprintf('Estimating MFCW distance errors based on linear model...'));

% 2FCW distance error at f1 vs. frequency difference
%     prepare matrix
linmodel.mfcw.derr_df = nan(settings.nfch, settings.np);
%     calculate
for i = 1 : settings.np
   linmodel.mfcw.derr_df(:, i) = pi/180*(arg.rho_d(ind.fch1, i) - arg.rho_d(:,i)) .* 3e8./(4*pi*(fch-fch(ind.fch1)));
end

% 2FCW distance error; fixed frequency difference vs. center frequency (f0, f0+df)
%     prepare cell array of matrices
for k = 1 : length(settings.mfcw.df12)
   linmodel.mfcw.derr_f0{k} = nan(settings.nfch, settings.np);
end
%     calculate
warning('off', 'MATLAB:interp1:NaNinY');
warning('off', 'MATLAB:chckxy:IgnoreNaN');
for i = 1 : settings.nfch
   for j = 1 : settings.np
      for k = 1 : length(settings.mfcw.df12)
         linmodel.mfcw.derr_f0{k}(i, j) = pi/180 * 3e8./(4*pi*settings.mfcw.df12(k)) * ...
            ( arg.rho_d(i, j) - interp1(fch, arg.rho_d(:, j), fch(i)+settings.mfcw.df12(k), 'linear','extrap') );
      end
   end
end
warning('on', 'MATLAB:interp1:NaNinY');
warning('on', 'MATLAB:chckxy:IgnoreNaN');

% 2FCW histogram of distance error over frequency
%     histogram setup
linmodel.mfcw.hist_derr = linspace(settings.mfcw.hist_minmax(1), settings.mfcw.hist_minmax(2), settings.mfcw.hist_n);
%     process distance error data
for j = 1 : settings.np
   for k = 1 : length(settings.mfcw.df12)
      % ... only for functional tags, P >= Pmin(unmod) 
      line_fch = linmodel.mfcw.derr_f0{k}(:, j);
      line_fch( pav(j) < pavopmin.all(:, end) ) = NaN;
      % create histogram and normalize
      linmodel.mfcw.hist_derr_f0{k}(j, :) = hist(line_fch, linmodel.mfcw.hist_derr);
      linmodel.mfcw.hist_derr_f0{k}(j, :) = linmodel.mfcw.hist_derr_f0{k}(j, :) / sum(linmodel.mfcw.hist_derr_f0{k}(j, :));
      % calculate expected value
      linmodel.mfcw.avg_derr_f0{k}(j) = nanmedian(line_fch);
   end
end


% *******************************************************************************************************
% plots
disp(sprintf('Creating plots...'));

% create templates for plot (group) names
% ... plot names do not have to be unique (will be numbered anyway)
plotnames = {};

% mute warnings (not all legend warnings have an ID => mute all)
oldwarn = warning('off', 'all');


% ind_f = interp1(fch, 1:settings.nfch, [850e6, 900e6, 950e6], 'nearest');
% rho_mod = squeeze(mag.rho(ind_f, 1, :));
% figure;
% subplot(2,1,1); hold on;
% plot(10*log10(pav)+30, squeeze(mag.rho(ind_f, 1, :)), '-');
% plot(10*log10(pav)+30, squeeze(mag.rho(ind_f, end, :)), '--');
% hold off; grid on; xlim([-15, 15]); ylim(xyzlimits(squeeze(mag.rho(ind_f, [1,settings.nm], :))));
% setlegend({'850MHz', '900MHz', '950MHz'}, 'SouthEast'); setlabels('', 'P_{av} [dBm]', '|S11|');
% subplot(2,1,2); hold on;
% plot(10*log10(pav)+30, squeeze(arg.rho(ind_f, 1, :)), '-');
% plot(10*log10(pav)+30, squeeze(arg.rho(ind_f, end, :)), '--');
% hold off; grid on; xlim([-15, 15]); ylim(xyzlimits(squeeze(arg.rho(ind_f, [1,settings.nm], :))));
% setlegend({'850MHz', '900MHz', '950MHz'}, 'NorthEast'); setlabels('', 'P_{av} [dBm]', 'arg(S11)');
% savefigure(gcf, 'tagchar_mod_600_0_preliminary', 'eps')
% 


% settings.figopt.visible = 'on'; return


% **************************************************
% general plots

% verification of interpolation and extrapolation
%     Ric(f, P)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_interp-ric']; %#ok<*AGROW>
surface(meas_tag_zic.fic/1e9, meas_tag_zic.pic, meas_tag_zic.ric , 'FaceColor','none');
surface(             fch/1e9,              pic,              ric', 'EdgeColor','none');
hold off; grid on; axis tight; view(-45, 30); set(gca, 'yScale','log', 'xMinorGrid', 'on');
setlabels('VERIFICATION: INTERPOLATION OF CHIP IMPEDANCE VS FREQUENCY/POWER', 'f [GHz]', 'P_{ic} [W]', 'R_{ic} [\Omega]');
setlegend({'original', 'modified'}, 'NorthEast');
%     Xic(f, P)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_interp-xic'];
surface(meas_tag_zic.fic/1e9, meas_tag_zic.pic, meas_tag_zic.xic , 'FaceColor','none');
surface(             fch/1e9,              pic,              xic', 'EdgeColor','none');
hold off; grid on; axis tight; view(-45, 30); set(gca, 'yScale','log', 'xMinorGrid','on');
setlabels('VERIFICATION: INTERPOLATION OF CHIP IMPEDANCE VS FREQUENCY/POWER', 'f [GHz]', 'P_{ic} [W]', 'X_{ic} [\Omega]');
setlegend({'original', 'modified'}, 'NorthEast');
%     Za(f)
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'general_interp-za'];
cplot(meas_tag_za.fa/1e9, meas_tag_za.ra, {'-',1,'blue'});
cplot(           fch/1e9,             ra, {}, {'o',3,'blue'});
cplot(meas_tag_za.fa/1e9, meas_tag_za.xa, {'-',1,'black'});
cplot(           fch/1e9,             xa, {}, {'o',3,'black'});
cplot(repmat(settings.fopmin+settings.foptol, 2, 1)/1e9, [0, 1.1*max([meas_tag_za.ra; meas_tag_za.xa])], {'-',1,'red'});
cplot(repmat(settings.fopmax-settings.foptol, 2, 1)/1e9, [0, 1.1*max([meas_tag_za.ra; meas_tag_za.xa])], {'-',1,'red'});
hold off; grid on; axis tight;
setlabels('VERIFICATION: INTERPOLATION OF ANTENNA IMPEDANCE VS FREQUENCY', 'f [GHz]', 'R_{a}, X_{a} [\Omega]');
setlegend({'R_{a} (original)', 'R_{a} (interpolated)',...
   'X_{a} (original)', 'X_{a} (interpolated)', 'area of operation'}, 'NorthWest');
%     |Zmod(P, mod)| in the middle of meas_tag_zmod.fic
meas_tag_zmod.zmod = complex(meas_tag_zmod.rmod, meas_tag_zmod.xmod);
ind_f = interp1(fch, 1:settings.nfch, meas_tag_zmod.fic(ceil(end/2)), 'nearest');
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_interp-rmod'];
plot3(meas_tag_zmod.pic, 100*ones(size(meas_tag_zmod.pic)), abs(meas_tag_zmod.zmod(:,(ceil(end/2))))', 'r', 'LineWidth', 5);
surface(pic, linspace(100,0,settings.nm), abs(squeeze(zmod(ind_f,:,:))), log10(abs(squeeze(zmod(ind_f,:,:)))), 'FaceColor','interp');
hold off; grid on; axis tight; view(-135, 50); set(gca, 'xScale','log', 'zScale','log');
setlabels('VERIFICATION: CHIP MODULATION IMPEDANCE VS POWER/MODULATION', 'P_{ic} [W]', 'modulation [%]',...
   sprintf('|Z_{mod}| [\\Omega] @ %.0fMHz', fch(ind_f)/1e6)); setlegend({'source data', 'extrapolated model'}, 'NorthEast');
%     Xmod(P, mod) in the middle of meas_tag_zmod.fic
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_interp-cmod'];
plot3(meas_tag_zmod.pic, 100*ones(size(meas_tag_zmod.pic)), 180/pi*angle(meas_tag_zmod.zmod(:,(ceil(end/2))))', 'r', 'LineWidth', 5);
surface(pic, linspace(100,0,settings.nm), 180/pi*angle(squeeze(zmod(ind_f,:,:))), 'FaceColor','interp');
hold off; grid on; axis tight; view(-135, 50); set(gca, 'xScale','log');
setlabels('VERIFICATION: CHIP MODULATION IMPEDANCE VS POWER/MODULATION', 'P_{ic} [W]', 'modulation [%]',...
   sprintf('arg(Z_{mod}) [deg] @ %.0fMHz', fch(ind_f)/1e6)); setlegend({'source data', 'extrapolated model'}, 'NorthEast');
%     Picmin(f)
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'general_interp-picmin'];
cplot(meas_tag_picmin.f/1e9, 10*log10(meas_tag_picmin.p)+30, {}, {'o',3,'blue'});
cplot(            fch/1e9, 10*log10(       picopmin)+30, {'-',1,'blue'});
hold off; grid on; xlim([min(fch),max(fch)]/1e9); ylim(xyzlimits(10*log10(meas_tag_picmin.p), 10*log10(picopmin))+30);
setlabels('VERIFICATION: INTERPOLATION OF MINIMUM OPERATIONAL CHIP POWER VS FREQUENCY', 'f [GHz]', 'P_{min} [dBm]');

% verification of rho_pav calculation (reverse check)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_check'];
%     real part
subplot(2,1,1); hold on;
cplot(err(:, 1)*100, {}, {'o',3,'black'});
cplot( ones(settings.checks,1)*settings.maxerr*100, {'-',1,'red'});
cplot( ones(settings.checks,1)*settings.maxavg*100, {'-',1,'green'});
cplot(-ones(settings.checks,1)*settings.maxerr*100, {'-',1,'red'});
cplot(-ones(settings.checks,1)*settings.maxavg*100, {'-',1,'green'});
hold off; grid on; axis tight; ylim(xyzlimits([settings.maxerr, -settings.maxerr]*100));
setlabels('RELATIVE MISMATCH Z_{ic} => \rho(f, Z_{mod}, P_{av}) => R_{ic}''',...
   'random check no.', 'relerr\{ R_{ic}, R_{ic}'' \} [%]');
setlegend({'error', 'error bound (for warning)', 'average error bound (for warning)'}, 'NorthEast');
%     imaginary part
subplot(2,1,2); hold on;
cplot(err(:, 2)*100, {}, {'o',3,'black'});
cplot( ones(settings.checks,1)*settings.maxerr*100, {'-',1,'red'});
cplot( ones(settings.checks,1)*settings.maxavg*100, {'-',1,'green'});
cplot(-ones(settings.checks,1)*settings.maxerr*100, {'-',1,'red'});
cplot(-ones(settings.checks,1)*settings.maxavg*100, {'-',1,'green'});
hold off; grid on; axis tight; ylim(xyzlimits([settings.maxerr, -settings.maxerr]*100));
setlabels('RELATIVE MISMATCH Z_{ic} => \rho(f, Z_{mod}, P_{av}) => X_{ic}''',...
   'random check no.', 'relerr\{ X_{ic}, X_{ic}'' \} [%]');
setlegend({'error', 'error bound (for warning)', 'average error bound (for warning)'}, 'NorthEast');

% position of largest errors
%     in magnitude plot
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'general_check-pos-mag'];
mesh(fch/1e9, pav, squeeze(mag.rho(:,   1, :))', 'EdgeColor','interp', 'FaceColor','none');
mesh(fch/1e9, pav, squeeze(mag.rho(:, end, :))', 'EdgeColor','interp', 'FaceColor','none');
plot3(fch(ind.err(1:settings.nwarperr,1))/1e9, pav(ind.err(1:settings.nwarperr,3)),...
   mag.rhoerr(:,1), 'k*', 'MarkerFaceColor', 'black', 'MarkerSize', 3);
plot3(fch(ind.err(1:settings.nwarperr,1))/1e9, pav(ind.err(1:settings.nwarperr,3)),...
   mag.rhoerr(:,2), 'ko', 'MarkerFaceColor', 'black', 'MarkerSize', 5);
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); axis tight; view([45, 30]);
setlabels(sprintf('POSITION OF %i LARGEST MISMATCHES Z_{ic} => \\rho(f, Z_{mod}, P_{av}) => Z_{ic}''', settings.nwarperr),...
   'f [GHz]', 'P_{av} [W]', '|\rho|');
setlegend({'modulated', 'unmodulated', 'small error (no problem)', 'large error'}, 'East');
%     in phase plot
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'general_check-pos-phase'];
mesh(fch/1e9, pav, squeeze(arg.rho(:, 1, :))',   'EdgeColor','interp', 'FaceColor','none');
mesh(fch/1e9, pav, squeeze(arg.rho(:, end, :))', 'EdgeColor','interp', 'FaceColor','none');
plot3(fch(ind.err(1:settings.nwarperr,1))/1e9, pav(ind.err(1:settings.nwarperr,3)),...
   arg.rhoerr(:,1), 'k*', 'MarkerFaceColor', 'black', 'MarkerSize', 3);
plot3(fch(ind.err(1:settings.nwarperr,1))/1e9, pav(ind.err(1:settings.nwarperr,3)),...
   arg.rhoerr(:,2), 'ko', 'MarkerFaceColor', 'black', 'MarkerSize', 5);
hold off; grid on; axis tight; set(gca, 'YScale','log', 'YTick',settings.pav_tick); view([45, 30]);
setlabels(sprintf('POSITION OF %i LARGEST MISMATCHES Z_{ic} => \\rho(f, Z_{mod}, P_{av}) => Z_{ic}''', settings.nwarperr),...
   'f [GHz]', 'P_{av} [W]', 'arg(\rho) [deg]');
setlegend({'modulated', 'unmodulated', 'small error (no problem)', 'large error'}, 'East');

% residual quadratic error (IIR bandpass optimization)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_iir-reserr'];
surface(linspace(100,0,settings.nm), pav, 20*log10(reserr)', 'EdgeColor','none', 'FaceColor','flat');
grid on; set(gca, 'YScale','log'); axis tight; xlim(xyzlimits([0,100])); ylim(xyzlimits(pav,'log'));
if ~isempty(settings.clim.reserr); set(gca, 'CLim',settings.clim.reserr); end; setcolorbar();
setlabels(sprintf('RESIDUAL QUADRATIC ERROR (COST FCN) OF IIR BANDPASS FITTING IN dB (MSE = %.1f dB)',...
   10*log10(nanmean(reserr(:).^2))), 'modulation [%]', 'P_{av} [W]', 'residual quadratic error [dB]');

% best IIR bandstop setup
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_iir-set'];
surface(linspace(100,0,settings.nm), pav, bestiirset', 'EdgeColor','none', 'FaceColor','flat');
grid on;  axis tight; setcolorbar(); xlim(xyzlimits([0,100])); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'CLim',[1,length(settings.iir)]);
 setlabels('BEST IIR BANDSTOP PARAMTER SET (IIR BANDPASS FITTING)',...
   'modulation [%]', 'P_{av} [W]', 'index for settings.iir cell array');

% needed time to fit IIR (IIR bandpass optimization)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_iir-time'];
surface(linspace(100,0,settings.nm), pav, t_fitiir', 'EdgeColor','none', 'FaceColor','flat');
grid on; set(gca, 'YScale','log'); xlim(xyzlimits([0,100])); ylim(xyzlimits(pav,'log'));
if ~isempty(settings.clim.t_fitiir); set(gca, 'CLim',settings.clim.t_fitiir); end; setcolorbar();
setlabels(sprintf('NEEDED TIME FOR IIR BANDPASS FITTING IN SECONDS (AVG = %.1f s)',...
   nanmean(t_fitiir(:))), 'modulation [%]', 'P_{av} [W]', 'fitting time [s]');

% influence of assembly tolerances on reflection coefficient magnitude
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'general_assembly'];
%     unmodulated
subplot(2,1,1);
surface(rat_at, cat_at*1e15, abs(squeeze(rho_at(1,:,:)))', 'EdgeColor','none');
grid on; axis tight; setcolorbar();
setlabels(sprintf('MAGNITUDE OF UNMODULATED REFLECTION COEFFICIENT @ %.0f MHz and %.1f dBm VS (DE)TUNING (PARALLEL)',...
   settings.fatm*1e-6, 10*log10(interp1(fch, picopmin, settings.fatm, 'nearest'))+30), 'R_{at} [\Omega]', 'C_{at} [fF]');
%     modulated
subplot(2,1,2);
surface(rat_at, cat_at*1e15, abs(squeeze(rho_at(2,:,:)))', 'EdgeColor','none');
grid on; axis tight; setcolorbar();
setlabels(sprintf('MAGNITUDE OF MODULATED REFLECTION COEFFICIENT @ %.0f MHz and %.1f dBm VS (DE)TUNING (SERIAL)',...
   settings.fatm*1e-6, 10*log10(interp1(fch, picopmin, settings.fatm, 'nearest'))+30), 'R_{at} [\Omega]', 'C_{at} [fF]');

% minimum available power necessary for tag operation
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'general_pavmin'];
cplot(fch/1e9, 10*log10(pavopmin.all(:,1))+30, {'-',1,'black'}); % for legend
cplot(fch/1e9, 10*log10(pavopmin.square50)+30, {'-',2,'red'});
cplot(fch/1e9, 10*log10(pavopmin.all)+30, {'-',1,'black'});
hold off; grid on; xlim([min(fch),max(fch)]/1e9); ylim(xyzlimits(10*log10(pavopmin.all)+30));
setlegend({'modulated through unmodulated', 'square mod. 50% duty cycle'}, 'NorthEast')
setlabels('MINIMUM AVAILABLE POWER NECESSARY FOR TAG OPERATION', 'f [GHz]', 'P_{av,min} [dBm]');


% **************************************************
% reflection coefficient non-extrapolated / extrapolated

% rho vs f an Pav for Zmod=min and Zmod=max (non-extrapolated)
%     magnitude
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_nonext-mag'];
surf(fch/1e9, pav, squeeze(mag.rho(:, 1, :))',   'EdgeColor','black', 'FaceColor','interp');
surf(fch/1e9, pav, squeeze(mag.rho(:, end, :))', 'EdgeColor','none',  'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick', settings.pav_tick); axis tight; view([20, 30]);
setlabels('MAGNITUDE OF REFLECTION COEFFICIENT (MODULATED / UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', '|\rho|');
setlegend({'modulated', 'unmodulated'}, 'NorthWest'); setcolorbar();
%     phase
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_nonext-phase'];
surf(fch/1e9, pav, squeeze(arg.rho(:, 1, :))',   'EdgeColor','black', 'FaceColor','interp');
surf(fch/1e9, pav, squeeze(arg.rho(:, end, :))', 'EdgeColor','none',  'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick', settings.pav_tick); axis tight; view([20, 30]);
setlabels('PHASE OF REFLECTION COEFFICIENT (MODULATED / UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'arg(\rho) [deg]');
setlegend({'modulated', 'unmodulated'}, 'NorthWest'); setcolorbar();

% rho vs f an Pav for Zmod=min and Zmod=max (extrapolated)
%     magnitude
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_ext-mag'];
surf(f/1e9, pav, squeeze(mag.rho_ext(:, 1, :))',   'EdgeColor','black', 'FaceColor','interp');
surf(f/1e9, pav, squeeze(mag.rho_ext(:, end, :))', 'EdgeColor','none',  'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'XTick',settings.fop_tick, 'YTick',settings.pav_tick); axis tight; xlim([settings.f_min, settings.f_max]); view([20, 30]);
setlabels('MAGNITUDE OF EXTRAPOLATED REFLECTION COEFFICIENT (MODULATED / UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', '|\rho|');
setlegend({'modulated', 'unmodulated'}, 'West'); setcolorbar();
%     phase
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_ext-phase'];
surf(f/1e9, pav, squeeze(arg.rho_ext(:, 1, :))',   'EdgeColor','black', 'FaceColor','interp');
surf(f/1e9, pav, squeeze(arg.rho_ext(:, end, :))', 'EdgeColor','none',  'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'XTick',settings.fop_tick, 'YTick', settings.pav_tick); axis tight; xlim([settings.f_min, settings.f_max]); view([20, 30]);
setlabels('PHASE OF EXTRAPOLATED REFLECTION COEFFICIENT (MODULATED / UNMODULATED)', 'f [GHz]', 'P_{av} [W]', 'arg(\rho) [deg]');
setlegend({'modulated', 'unmodulated'}, 'NorthEast'); setcolorbar();

% comparison extrapolated/nonextrapolated vs f an Pav for Zmod=max (unmodulated)
%     magnitude
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_comp'];
surf(  f/1e9, pav, squeeze(mag.rho_ext(:, end, :))', 'EdgeColor','black', 'FaceColor','none');
surf(fch/1e9, pav, squeeze(mag.rho    (:, end, :))', 'EdgeColor','none',  'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick', settings.pav_tick); axis tight; xlim([fch(1)/1e9-0.01, fch(end)/1e9+0.01]); view([20, 50]);
setlabels('MAGNITUDE OF UNMODULATED REFLECTION COEFFICIENT', 'f [GHz]', 'P_{av} [W]', '|\rho|');
setlegend({'extrapolated', 'non-extrapolated'}, 'NorthWest'); setcolorbar();
%     SPECIAL MODIFICATION FOR PRESENTATIONS AND THESIS
lp = length(pav); ind_p = [1:5,6:3:lp, lp];
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_comp'];
surf(  f/1e9, pav(ind_p), squeeze(mag.rho_ext(:, end, ind_p))', 'EdgeColor','black', 'FaceColor','none');
surf(fch/1e9, pav, squeeze(mag.rho    (:, end, :))', 'EdgeColor','none',  'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick', settings.pav_tick); axis tight; xlim([fch(1)/1e9-0.01, fch(end)/1e9+0.01]); view([20, 50]);
setlabels('', 'f [GHz]', 'P_{av} [W]', '|\rho|'); setlegend({'extrapolated', 'non-extrapolated'}, 'NorthWest'); setcolorbar();

% magnitude of rho vs f an Pav (extrapolated, xy-plane)
%     for Zmod=max (unmodulated)
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_ext-xy-mod'];
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1); 
surface(f/1e9, pav, squeeze(mag.rho_ext(:, end, :))', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'XTick',settings.fop_tick, 'YTick',settings.pav_tick); axis tight; xlim([settings.f_min, settings.f_max]);
setlabels('MAGNITUDE OF EXTRAPOLATED REFLECTION COEFFICIENT (UNMODULATED)', 'f [GHz]', 'P_{av} [W]', '|\rho|');
setlegend({'P_{min}'}, 'NorthWest'); setcolorbar();
%     for Zmod=min (modulated)
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_ext-xy-unmod'];
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.np, settings.nfch), 'Color','k', 'LineWidth',1);
surface(f/1e9, pav, squeeze(mag.rho_ext(:, 1, :))', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'XTick',settings.fop_tick, 'YTick',settings.pav_tick); axis tight; xlim([settings.f_min, settings.f_max]);
setlabels('MAGNITUDE OF EXTRAPOLATED REFLECTION COEFFICIENT (MODULATED)', 'f [GHz]', 'P_{av} [W]', '|\rho|');
setlegend({'P_{min}'}, 'NorthWest'); setcolorbar();

% rho vs |Zmod| and P (non-extrapolated)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'rho_zmod'];
%     magnitude
subplot(2,1,1);
ind.temp = findzeros([1, all(isnan(squeeze(mag.rho(ind.f1, :, :))), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
plot3(linspace(100,0,settings.nm), squeeze(pavopmin.all(ind.f1,:)), 1.1*ones(settings.nm,1), 'Color','k', 'LineWidth',1);
surface(linspace(100,0,settings.nm), pav, squeeze(mag.rho(ind.f1, :, :))', 'EdgeColor','none', 'FaceColor','interp');
grid on; set(gca, 'YScale','log', 'YTick', settings.pav_tick); axis tight; 
xlim(xyzlimits([0,100])); ylim(xyzlimits([pav(ind.temp(1)), pav(ind.temp(2))], 'log')); view(0,90);
setlabels(sprintf('MAGNITUDE OF NON-EXTRAPOLATED REFLECTION COEFFICIENT @ %.0f MHz',  f(ind.f1)/1e6),...
   'modulation [%]', 'P_{av} [W]', '|\rho|');
setlegend({'P_{min}'}, 'SouthWest'); setcolorbar();
%     phase
subplot(2,1,2);
ind.temp = findzeros([1, all(isnan(squeeze(arg.rho(ind.f1, :, :))), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
plot3(linspace(100,0,settings.nm), squeeze(pavopmin.all(ind.f1,:)),...
   1.1*max(max(unwrap(squeeze(arg.rho(ind.f1, :, :)), 180, 1)))*ones(settings.nm,1), 'Color','k', 'LineWidth',1);
surface(linspace(100,0,settings.nm), pav, unwrap(squeeze(arg.rho(ind.f1, :, :)), 180, 1)', 'EdgeColor','none', 'FaceColor','interp');
grid on; set(gca, 'YScale','log', 'YTick', settings.pav_tick); axis tight; 
xlim(xyzlimits([0,100])); ylim(xyzlimits([pav(ind.temp(1)), pav(ind.temp(2))], 'log')); view(0,90);
setlabels(sprintf('PHASE OF NON-EXTRAPOLATED REFLECTION COEFFICIENT @ %.0f MHz IN DEGREE',  f(ind.f1)/1e6),...
   'modulation [%]', 'P_{av} [W]', 'arg(\rho)');
setlegend({'P_{min}'}, 'SouthWest'); setcolorbar();

% delta_rho (maximum stroke) vs f an Pav
%     magnitude
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_stroke-mag'];
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1); 
surface(fch/1e9, pav, squeeze(mag.rho(:, 1, :))' - squeeze(mag.rho(:, end, :))', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log'); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
setlabels('DIFFERENCE IN REFLECTION COEFFICIENT MAGNITUDE (MODULATED / UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', '|\rho|');
%     magnitude with fixed limits (for comparisons)
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_stroke-mag-comp'];
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1); 
surface(fch/1e9, pav, squeeze(mag.rho(:, 1, :))' - squeeze(mag.rho(:, end, :))', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'CLim',settings.clim.delta_rho); 
xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
setlabels('DIFFERENCE IN REFLECTION COEFFICIENT MAGNITUDE (MODULATED / UNMODULATED) FOR COMPARISON',...
   'f [GHz]', 'P_{av} [W]', '|\rho|');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
%     phase
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_stroke-phase'];
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1); 
surface(fch/1e9, pav, squeeze(arg.rho(:, 1, :))' - squeeze(arg.rho(:, end, :))', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log'); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
setlabels('DIFFERENCE IN REFLECTION COEFFICIENT PHASE (MODULATED / UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'arg(\rho)');

% delta_rho (maximum stroke) vs f an Pav (extrapolated; only magnitude important)
figure('Visible', settings.figopt.visible); hold on; plotnames = [plotnames, 'rho_ext-stroke-mag'];
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1); 
surface(f/1e9, pav, squeeze(mag.rho_ext(:, 1, :))' - squeeze(mag.rho_ext(:, end, :))', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'XTick',settings.fop_tick); axis tight; xlim([settings.f_min, settings.f_max]);
setlabels('DIFFERENCE IN EXTRAP. REFLECTION COEFFICIENT MAGNITUDE (MODULATED / UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', '|\rho|');
 setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();


% **************************************************
% linear model

% linear model: center rho_c vs f an Pav
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_center'];
%     magnitude
subplot(2,1,1); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, mag.rho_c', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('LINEAR MODEL: CENTER VALUE \rho_c (MAGNITUDE)', 'f [GHz]', 'P_{av} [W]', '|\rho_c|');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
%     phase
subplot(2,1,2); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, arg.rho_c', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('LINEAR MODEL: CENTER VALUE \rho_c (PHASE IN DEGREE)', 'f [GHz]', 'P_{av} [W]', 'arg(\rho_c) [deg]');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();

% linear model: difference rho_d vs f an Pav
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_diff'];
%     magnitude
subplot(2,1,1); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, mag.rho_d', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('LINEAR MODEL: DIFFERENCE \rho_{\Delta} (MAGNITUDE)', 'f [GHz]', 'P_{av} [W]', '|\rho_{\Delta}|');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
%     phase
subplot(2,1,2); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, arg.rho_d', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('LINEAR MODEL: DIFFERENCE \rho_{\Delta} (PHASE IN DEGREE)', 'f [GHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [deg]');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();

% linear model:center rho_c vs f an Pav with fixed color limits for comparison
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_center-comp'];
%     magnitude
subplot(2,1,1); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, mag.rho_c', 'EdgeColor','none', 'FaceColor','interp');
grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('LINEAR MODEL: CENTER VALUE \rho_c (MAGNITUDE) FOR COMPARISON', 'f [GHz]', 'P_{av} [W]', '|\rho_c|');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
if ~isempty(settings.clim.rho_c_mag); set(gca, 'CLim', settings.clim.rho_c_mag); end;
%     phase
subplot(2,1,2); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, arg.rho_c', 'EdgeColor','none', 'FaceColor','interp');
grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('LINEAR MODEL: CENTER VALUE \rho_c (PHASE IN DEGREE) FOR COMPARISON', 'f [GHz]', 'P_{av} [W]', 'arg(\rho_c) [deg]');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
if ~isempty(settings.clim.rho_c_arg); set(gca, 'CLim', settings.clim.rho_c_arg); end;

% linear model: difference rho_d vs f an Pav with fixed color limits for comparison
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_diff-comp'];
%     magnitude
subplot(2,1,1); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, mag.rho_d', 'EdgeColor','none', 'FaceColor','interp');
grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('LINEAR MODEL: DIFFERENCE \rho_{\Delta} (MAGNITUDE) FOR COMPARISON', 'f [GHz]', 'P_{av} [W]', '|\rho_{\Delta}|');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
if ~isempty(settings.clim.rho_d_mag); set(gca, 'CLim', settings.clim.rho_d_mag); end;
%     phase
subplot(2,1,2); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, arg.rho_d', 'EdgeColor','none', 'FaceColor','interp');
grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('LINEAR MODEL: DIFFERENCE \rho_{\Delta} (PHASE IN DEGREE) FOR COMPARISON', 'f [GHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [deg]');
setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
if ~isempty(settings.clim.rho_d_arg); set(gca, 'CLim', settings.clim.rho_d_arg); end;

% linear model: drho_d/df vs f an Pav
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_grad'];
%        magnitude
subplot(2,1,1); hold on;
surface(fch_drho/1e9, pav, mag.drho_d'*1e6, 'EdgeColor','none', 'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(mag.drho_d(:))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim', [min(min(mag.drho_d)), max(max(mag.drho_d))]*1e6);
setlabels('LINEAR MODEL: GRADIENT OF \rho_{\Delta} MAGNITUDE [1/MHz]', 'f [GHz]', 'P_{av} [W]', 'd/df arg(\rho) [1/MHz]');
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
%        phase
subplot(2,1,2); hold on;
surface(fch_drho/1e9, pav, arg.drho_d'*1e6, 'EdgeColor','none',  'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(arg.drho_d(:))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim', [min(min(arg.drho_d)), max(max(arg.drho_d))]*1e6);
setlabels('LINEAR MODEL: GRADIENT OF \rho_{\Delta} PHASE [deg/MHz]', 'f [GHz]', 'P_{av} [W]', 'd/df arg(\rho) [deg/MHz]');
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');

% linear model: drho_d/df vs f an Pav using zoomed colormap for better visibility of smaller changes
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_grad-zoom'];
%        magnitude
subplot(2,1,1); hold on;
surface(fch_drho/1e9, pav, mag.drho_d'*1e6, 'EdgeColor','none', 'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(mag.drho_d(:))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim',...
   [min(min(mag.drho_d(ind.fopmin:ind.fopmax,:))), max(max(mag.drho_d(ind.fopmin:ind.fopmax,:)))]*1e6);
 xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('[zoomed color limits] LINEAR MODEL: GRADIENT OF \rho_{\Delta} MAGNITUDE [1/MHz]', 'f [GHz]', 'P_{av} [W]', 'd/df arg(\rho) [1/MHz]');
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
%        phase
subplot(2,1,2); hold on;
surface(fch_drho/1e9, pav, arg.drho_d'*1e6, 'EdgeColor','none',  'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(arg.drho_d(:))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim',...
   [min(min(arg.drho_d(ind.fopmin:ind.fopmax,:))), max(max(arg.drho_d(ind.fopmin:ind.fopmax,:)))]*1e6);
 xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('[zoomed color limits] LINEAR MODEL: GRADIENT OF \rho_{\Delta} PHASE [deg/MHz]', 'f [GHz]', 'P_{av} [W]', 'd/df arg(\rho) [deg/MHz]');
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');

% linear model: separate plot for drho_d (phase)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_diff-phase']; hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, arg.rho_d', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('', 'f_c [GHz]', 'P_{av} [W]'); setlegend({'P_{av,min} (mod, unmod)'}, 'NorthEast'); setcolorbar();
disp('warning: modified caption');
% linear model: separate plot for drho_d/df (phase)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_grad-phase']; hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch_drho/1e9, pav, arg.drho_d'*1e6, 'EdgeColor','none',  'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim',...
   [min(min(arg.drho_d(ind.fopmin:ind.fopmax,:))), max(max(arg.drho_d(ind.fopmin:ind.fopmax,:)))]*1e6);
 xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels('', 'f_c [GHz]', 'P_{av} [W]'); setlegend({'P_{av,min} (mod, unmod)'}, 'NorthEast'); setcolorbar();



% **************************************************
% MFCW (multi-frequency continuous-wave ranging)

% MFCW: distance error for carrier at fch, one secondary carrier at fch+df12 (2FCW)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'mfcw_derr-f0'];
%     at first frequency difference
subplot(2,1,1); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 101*max(max(linmodel.mfcw.derr_f0{1}))*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 101*max(max(linmodel.mfcw.derr_f0{1}))*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, 100*linmodel.mfcw.derr_f0{1}', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels(sprintf('LINEAR MODEL: DISTANCE ERROR [cm] FOR 2FCW, FREQ. DIFFERENCE f_1-f_0 = %.1f MHz', settings.mfcw.df12(1)/1e6),...
   'f_0 [GHz] (main carrier)', 'P_{av} [W]', 'd_{err} [cm]');
setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
if ~isempty(settings.clim.mfcw_derr1); set(gca, 'CLim', 100*settings.clim.mfcw_derr1); end;
%     at second frequency difference
subplot(2,1,2); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 101*max(max(linmodel.mfcw.derr_f0{2}))*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 101*max(max(linmodel.mfcw.derr_f0{2}))*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, 100*linmodel.mfcw.derr_f0{2}', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels(sprintf('LINEAR MODEL: DISTANCE ERROR [cm] FOR 2FCW, FREQ. DIFFERENCE f_1-f_0 = %.1f MHz', settings.mfcw.df12(2)/1e6),...
   'f_0 [GHz] (main carrier)', 'P_{av} [W]', 'd_{err} [cm]');
setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
if ~isempty(settings.clim.mfcw_derr1); set(gca, 'CLim', 100*settings.clim.mfcw_derr1); end;

% MFCW: histogram for distance error; main carrier at fch, one secondary carrier at fch+df12 (2FCW)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'mfcw_hist-derr-f0'];
%     at first frequency difference
subplot(2,1,1); hold on;
plot3(pav, 100*linmodel.mfcw.avg_derr_f0{1}, 101*max(max(linmodel.mfcw.hist_derr_f0{1}))*ones(settings.np,1), 'Color','k', 'LineWidth',1);
plot3(pav, zeros(size(pav)), 101*max(max(linmodel.mfcw.hist_derr_f0{1}))*ones(settings.np,1), ':', 'Color','k', 'LineWidth',1);
surface(pav, 100*linmodel.mfcw.hist_derr, 100*linmodel.mfcw.hist_derr_f0{1}', 'FaceColor','interp', 'EdgeColor','none');
hold off; grid on; set(gca, 'xScale','log', 'xTick',settings.pav_tick);
xlim(xyzlimits(pav,'log')); ylim(xyzlimits(100*linmodel.mfcw.hist_derr));
setlabels(sprintf('LINEAR MODEL: DISTANCE ERROR PDF [%%] FOR 2FCW, FREQ. DIFFERENCE f_1-f_0 = %.1f MHz', settings.mfcw.df12(1)/1e6),...
   'P_{av} [W]', 'd_{err} [cm]', 'PDF(d_{err}) @ P_{av} [%%]');
setlegend({'median'}, 'SouthWest'); setcolorbar();
if ~isempty(settings.clim.mfcw_derr_hist); set(gca, 'CLim', settings.clim.mfcw_derr_hist); end;
%     at second frequency difference
subplot(2,1,2); hold on;
plot3(pav, 100*linmodel.mfcw.avg_derr_f0{2}, 101*max(max(linmodel.mfcw.hist_derr_f0{2}))*ones(settings.np,1), 'Color','k', 'LineWidth',1);
plot3(pav, zeros(size(pav)), 101*max(max(linmodel.mfcw.hist_derr_f0{2}))*ones(settings.np,1), ':', 'Color','k', 'LineWidth',1);
surface(pav, 100*linmodel.mfcw.hist_derr, 100*linmodel.mfcw.hist_derr_f0{2}', 'FaceColor','interp', 'EdgeColor','none');
hold off; grid on; set(gca, 'xScale','log', 'xTick',settings.pav_tick);
xlim(xyzlimits(pav,'log')); ylim(xyzlimits(100*linmodel.mfcw.hist_derr));
setlabels(sprintf('LINEAR MODEL: DISTANCE ERROR PDF [%%] FOR 2FCW, FREQ. DIFFERENCE f_1-f_0 = %.1f MHz', settings.mfcw.df12(2)/1e6),...
   'P_{av} [W]', 'd_{err} [cm]', 'PDF(d_{err}) @ P_{av} [%%]');
setlegend({'median'}, 'SouthWest'); setcolorbar();
if ~isempty(settings.clim.mfcw_derr_hist); set(gca, 'CLim', settings.clim.mfcw_derr_hist); end;

% MFCW: distance error for carrier at settings.f1, one secondary carrier (2FCW)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'mfcw_derr-df'];
%     full span colormap
subplot(2,1,1); hold on;
plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*max(max(linmodel.mfcw.derr_df))*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*max(max(linmodel.mfcw.derr_df))*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e9, pav, linmodel.mfcw.derr_df', 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
setlabels(sprintf('LINEAR MODEL: DISTANCE ERROR [m] FOR 2FCW, MAIN CARRIER AT f_0 = %.1f MHz (f: SECONDARY CARRIER)', fch(ind.fch1)/1e6),...
   'f [GHz]', 'P_{av} [W]', 'd_{err} [m]');
setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
if ~isempty(settings.clim.mfcw_derr1); set(gca, 'CLim', settings.clim.mfcw_derr1); end; setcolorbar();
%     zoomed for better visibility
ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f1-settings.fd, 'nearest', 'extrap'));
ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f1+settings.fd, 'nearest', 'extrap'));
ind.temp = findzeros([1, all(isnan(linmodel.mfcw.derr_df(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
subplot(2,1,2); hold on;
plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 110*max(max(linmodel.mfcw.derr_df))*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 110*max(max(linmodel.mfcw.derr_df))*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
surface(fch/1e6, pav, linmodel.mfcw.derr_df'*100, 'EdgeColor','none', 'FaceColor','interp');
hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick);
xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim(xyzlimits([pav(ind.temp(1)), pav(ind.temp(2))],'log'));
setlabels(sprintf('[zoomed] LINEAR MODEL: DISTANCE ERROR [cm] FOR 2FCW, MAIN CARRIER AT f_0 = %.1f MHz (f: SECONDARY CARRIER)', settings.f1/1e6),...
   'f [MHz]', 'P_{av} [W]', 'd_{err} [cm]');
setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
if ~isempty(settings.clim.mfcw_derr2)
   set(gca, 'CLim', settings.clim.mfcw_derr2);
else
   set(gca, 'CLim', [min(min(linmodel.mfcw.derr_df(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.mfcw.derr_df(ind.fch_d1:ind.fch_d2,:)))]*100);
end


% **************************************************
% gradients (except for linear model)

% magnitude of drho/df vs f an Pav for Zmod=min and Zmod=max (non-extrapolated)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'grad_rho-mag'];
%        Zmod=min (modulated)
subplot(2,1,1); hold on;
surface(fch_drho/1e9, pav, squeeze(mag.drho(:, 1, :))'*1e6, 'EdgeColor','none', 'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(min(mag.drho(:, 1, :)))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; setcolorbar(); xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim', [min(min(mag.drho(:, 1, :))), max(max(mag.drho(:, 1, :)))]*1e6);
setlegend({'gradient', 'background: support'}, 'SouthWest');
setlabels('GRADIENT OF REFLECTION COEFFICIENT MAGNITUDE [1/MHz] (MODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'd/df |rho| [1/MHz]');
%        Zmod=max (unmodulated)
subplot(2,1,2); hold on;
surface(fch_drho/1e9, pav, squeeze(mag.drho(:, end, :))'*1e6, 'EdgeColor','none',  'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(min(mag.drho(:, end, :)))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim', [min(min(mag.drho(:, end, :))), max(max(mag.drho(:, end, :)))]*1e6);
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
setlabels('GRADIENT OF REFLECTION COEFFICIENT MAGNITUDE [1/MHz] (UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'd/df |rho| [1/MHz]');

% magnitude of drho/df vs f an Pav for Zmod=min and Zmod=max (non-extrapolated)
% using colormap of drhoe/df for better visibility of smaller changes
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'grad_rho-mag-zoom'];
%        Zmod=min (modulated)
subplot(2,1,1); hold on;
surface(fch_drho/1e9, pav, squeeze(mag.drho(:, 1, :))'*1e6, 'EdgeColor','none', 'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(min(mag.drho(:, 1, :)))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,1)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim', [min(min(mag.drhoe(:, 1, :))), max(max(mag.drhoe(:, 1, :)))]*1e6);
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
setlabels('[zoomed color limits] GRADIENT OF REFLECTION COEFFICIENT MAGNITUDE [1/MHz] (MODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'd/df |rho| [1/MHz]');
%        Zmod=max (unmodulated)
subplot(2,1,2); hold on;
surface(fch_drho/1e9, pav, squeeze(mag.drho(:, end, :))'*1e6, 'EdgeColor','none',  'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(min(mag.drho(:, end, :)))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim',...
   [min(min(mag.drhoe(ind.fopmin:ind.fopmax, end, :))), max(max(mag.drhoe(ind.fopmin:ind.fopmax, end, :)))]*1e6);
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
setlabels('[zoomed color limits] GRADIENT OF REFLECTION COEFFICIENT MAGNITUDE [1/MHz] (UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'd/df |rho| [1/MHz]');

% phase of drho/df vs f an Pav for Zmod=min and Zmod=max (non-extrapolated)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'grad_rho-phase'];
%        Zmod=min (modulated)
subplot(2,1,1); hold on;
surface(fch_drho/1e9, pav, squeeze(arg.drho(:, 1, :))'*1e6, 'EdgeColor','none', 'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(min(arg.drho(:, 1, :)))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,1)), 100*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim', [min(min(arg.drho(:, 1, :))), max(max(arg.drho(:, 1, :)))]*1e6);
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
setlabels('GRADIENT OF REFLECTION COEFFICIENT PHASE [deg/MHz] (MODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'd/df arg(\rho) [deg/MHz]');
%        Zmod=max (unmodulated)
subplot(2,1,2); hold on;
surface(fch_drho/1e9, pav, squeeze(arg.drho(:, end, :))'*1e6, 'EdgeColor','none',  'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(min(arg.drho(:, end, :)))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 100*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim', [min(min(arg.drho(:, end, :))), max(max(arg.drho(:, end, :)))]*1e6);
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
setlabels('GRADIENT OF REFLECTION COEFFICIENT PHASE [deg/MHz] (UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'd/df arg(\rho) [deg/MHz]');

% phase of drho/df vs f an Pav for Zmod=min and Zmod=max (non-extrapolated)
% using colormap of drhoe/df for better visibility of smaller changes
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'grad_rho-phase-zoom'];
%        Zmod=min (modulated)
subplot(2,1,1); hold on;
surface(fch_drho/1e9, pav, squeeze(arg.drho(:, 1, :))'*1e6, 'EdgeColor','none', 'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(min(arg.drho(:, 1, :)))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,1)), 100*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim',...
   [min(min(arg.drhoe(ind.fopmin:ind.fopmax, 1, :))), max(max(arg.drhoe(ind.fopmin:ind.fopmax, 1, :)))]*1e6);
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
setlabels('[zoomed color limits] GRADIENT OF REFLECTION COEFFICIENT PHASE [deg/MHz] (MODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'd/df arg(\rho) [deg/MHz]');
%        Zmod=max (unmodulated)
subplot(2,1,2); hold on;
surface(fch_drho/1e9, pav, squeeze(arg.drho(:, end, :))'*1e6, 'EdgeColor','none',  'FaceColor','interp');
surface(fch     /1e9, pav, background.f_pav_nan'-1+min(min(arg.drho(:, end, :)))*1e6,...
   repmat(-background.f_pav', [1,1,3]), 'EdgeColor','interp', 'FaceColor','none');
plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 100*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim(xyzlimits(pav,'log'));
set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'Clim',...
   [min(min(arg.drhoe(ind.fopmin:ind.fopmax, end, :))), max(max(arg.drhoe(ind.fopmin:ind.fopmax, end, :)))]*1e6);
setcolorbar(); setlegend({'gradient', 'background: support'}, 'SouthWest');
setlabels('[zoomed color limits] GRADIENT OF REFLECTION COEFFICIENT PHASE [deg/MHz] (UNMODULATED)',...
   'f [GHz]', 'P_{av} [W]', 'd/df arg(\rho) [deg/MHz]');


% **************************************************
% example filter characteristics

% different powers mod/unmod
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'filter_mod'];
%     magnitude
subplot(2,1,1); hold on;
cplot(fch/1e9, squeeze(mag.rho(:,1,  ind.pav1)), {'--',1,'black'});
cplot(fch/1e9, squeeze(mag.rho(:,end,ind.pav1)), {'-', 1,'black'});
cplot(fch/1e9, squeeze(mag.rho(:,1,  ind.pav2)), {'--',1,'blue' });
cplot(fch/1e9, squeeze(mag.rho(:,end,ind.pav2)), {'-', 1,'blue' });
hold off; grid on; xlim([settings.fchmin, settings.fchmax]/1e9); ylim(xyzlimits(mag.rho(:,[1;end],[ind.pav1;ind.pav2])));
setlabels('EXAMPLE REFLECTION FILTER FREQUENCY RESPONSE: MAGNITUDE', 'f [GHz]', '|\rho|');
setlegend({...
   sprintf('P_{av}=%.1e W, mod', pav(ind.pav1)), sprintf('P_{av}=%.1e W, unmod', pav(ind.pav1)),...
   sprintf('P_{av}=%.1e W, mod', pav(ind.pav2)), sprintf('P_{av}=%.1e W, unmod', pav(ind.pav2))},...
   'SouthEast');
%     phase
subplot(2,1,2); hold on;
cplot(fch/1e9, squeeze(arg.rho(:,1,  ind.pav1)), {'--',1,'black'});
cplot(fch/1e9, squeeze(arg.rho(:,end,ind.pav1)), {'-', 1,'black'});
cplot(fch/1e9, squeeze(arg.rho(:,1,  ind.pav2)), {'--',1,'blue' });
cplot(fch/1e9, squeeze(arg.rho(:,end,ind.pav2)), {'-', 1,'blue' });
hold off; grid on; xlim([settings.fchmin, settings.fchmax]/1e9); ylim(xyzlimits(arg.rho(:,[1;end],[ind.pav1;ind.pav2])));
setlabels('EXAMPLE REFLECTION FILTER FREQUENCY RESPONSE: PHASE', 'f [GHz]', ' arg(\rho) [deg]');
setlegend({...
   sprintf('P_{av}=%.1e W, mod', pav(ind.pav1)), sprintf('P_{av}=%.1e W, unmod', pav(ind.pav1)),...
   sprintf('P_{av}=%.1e W, mod', pav(ind.pav2)), sprintf('P_{av}=%.1e W, unmod', pav(ind.pav2))},...
   'NorthEast');

% extrapolated vs. non-extrapolated (unmod)
figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'filter_extrap'];
%     magnitude
subplot(2,1,1); hold on;
cplot(  f/1e9, squeeze(mag.rho_ext(:,end,ind.pav1)), {'-',1,'black'});
cplot(fch/1e9, squeeze(mag.rho    (:,end,ind.pav1)),              {}, {'o',3,'black'});
cplot(  f/1e9, squeeze(mag.rho_ext(:,end,ind.pav2)), {'-',1, 'blue'});
cplot(fch/1e9, squeeze(mag.rho    (:,end,ind.pav2)),              {}, {'o',3, 'blue'});
hold off; grid on; xlim([settings.f_min, settings.f_max]); set(gca, 'XTick',settings.fop_tick);
ylim(xyzlimits([mag.rho_ext(:,end,ind.pav1); mag.rho_ext(:,end,ind.pav2)])); 
setlabels('REFLECTION FILTER FREQUENCY RESPONSE EXTRAP. VS NONEXTRAP.: MAGNITUDE (UNMOD)', 'f [GHz]', '|\rho|');
setlegend({...
   sprintf('extrapolated (P_{av}=%.1e W)', pav(ind.pav1)), 'non-extrapolated ',...
   sprintf('extrapolated (P_{av}=%.1e W)', pav(ind.pav2)), 'non-extrapolated'}, 'SouthEast');
%     phase
subplot(2,1,2); hold on;
cplot(  f/1e9, squeeze(arg.rho_ext(:,end,ind.pav1)), {'-',1,'black'});
cplot(fch/1e9, squeeze(arg.rho    (:,end,ind.pav1)),              {}, {'o',3,'black'});
cplot(  f/1e9, squeeze(arg.rho_ext(:,end,ind.pav2)), {'-',1, 'blue'});
cplot(fch/1e9, squeeze(arg.rho    (:,end,ind.pav2)),              {}, {'o',3, 'blue'});
hold off; grid on; xlim([settings.f_min, settings.f_max]); set(gca, 'XTick',settings.fop_tick);
ylim(xyzlimits([arg.rho_ext(:,end,ind.pav1); arg.rho_ext(:,end,ind.pav2)]));
setlabels('REFLECTION FILTER FREQUENCY RESPONSE EXTRAP. VS NONEXTRAP.: PHASE (UNMOD)', 'f [GHz]', ' arg(\rho) [deg]');
setlegend({...
   sprintf('extrapolated (P_{av}=%.1e W)', pav(ind.pav1)), 'non-extrapolated ',...
   sprintf('extrapolated (P_{av}=%.1e W)', pav(ind.pav2)), 'non-extrapolated'}, 'NorthEast');


% **************************************************
% reactivate deactivated warnings
warning(oldwarn);


% *******************************************************************************************************
% save plots to file
disp(sprintf('Saving plots...'));

% check if all plots have been named
if gcf ~= length(plotnames)
   error('Mismatch between number of figures and length of plotnames array.');
end

% save and close figures
%   create new directory, meanwhile deactivate warning "Removed '/blah' from the MATLAB path..."
warning('off', 'MATLAB:RMDIR:RemovedFromPath');
[status, message, messageid] = rmdir(fullfile(settings.plotfolder, sprintf('tagchar_mod%s', settings.suffix)), 's'); %#ok<NASGU>
[status, message, messageid] = mkdir(fullfile(settings.plotfolder, sprintf('tagchar_mod%s', settings.suffix)));
warning('on', 'MATLAB:RMDIR:RemovedFromPath');
%   save figures
for i = 1 : gcf
   if isempty(plotnames{i})
      savefigure(i, fullfile(settings.plotfolder, sprintf('tagchar_mod%s', settings.suffix),...
         sprintf('tagchar_mod%s_%02i', settings.suffix, i)), settings.figopt.type, settings.figopt);
   else
      savefigure(i, fullfile(settings.plotfolder, sprintf('tagchar_mod%s', settings.suffix),...
         sprintf('tagchar_mod%s_%02i_%s', settings.suffix, i, plotnames{i})),...
         settings.figopt.type, settings.figopt);
   end
   close(i);
end




% *******************************************************************************************************
% *******************************************************************************************************
% OLD AND RUSTY

% % linear model:  phase of center rho_c vs f an Pav zoomed around some frequency values [1/2]
% figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear'];
% %     around settings.f1
% subplot(2,1,1); hold on;
% ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f1-settings.fd, 'nearest', 'extrap'));
% ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f1+settings.fd, 'nearest', 'extrap'));
% ind.temp = findzeros([1, all(isnan(arg.rho_c(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
% plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% surface(fch/1e6, pav, linmodel.darg.rho_c1', 'EdgeColor','none', 'FaceColor','interp');
% hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick,...
%    'CLim', [min(min(linmodel.darg.rho_c1(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.darg.rho_c1(ind.fch_d1:ind.fch_d2,:)))]);
% xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim([pav(ind.temp(1)), pav(ind.temp(2))]);
% setlabels(sprintf('[zoomed] LINEAR MODEL: PHASE DIFFERENCE OF \\rho_c TO f = %.1f MHz (IN DEGREE)', fch(ind.fch1)/1e6),...
%    'f [MHz]', 'P_{av} [W]', 'arg(\rho_c) [deg]');
% setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
% %     around settings.f2
% subplot(2,1,2); hold on;
% ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f2-settings.fd, 'nearest', 'extrap'));
% ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f2+settings.fd, 'nearest', 'extrap'));
% ind.temp = findzeros([1, all(isnan(arg.rho_c(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
% plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% surface(fch/1e6, pav, linmodel.darg.rho_c2', 'EdgeColor','none', 'FaceColor','interp');
% hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick,...
%    'CLim', [min(min(linmodel.darg.rho_c2(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.darg.rho_c2(ind.fch_d1:ind.fch_d2,:)))]);
% xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim([pav(ind.temp(1)), pav(ind.temp(2))]);
% setlabels(sprintf('[zoomed] LINEAR MODEL: PHASE DIFFERENCE OF \\rho_c TO f = %.1f MHz (IN DEGREE)', fch(ind.fch2)/1e6),...
%    'f [MHz]', 'P_{av} [W]', 'arg(\rho_c) [deg]');
% setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
% 
% % linear model:  phase of center rho_c vs f an Pav zoomed around some frequency values [2/2]
% figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear'];
% %     around settings.f3
% subplot(2,1,1); hold on;
% ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f3-settings.fd, 'nearest', 'extrap'));
% ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f3+settings.fd, 'nearest', 'extrap'));
% ind.temp = findzeros([1, all(isnan(arg.rho_c(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
% plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% surface(fch/1e6, pav, linmodel.darg.rho_c3', 'EdgeColor','none', 'FaceColor','interp');
% hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick,...
%    'CLim', [min(min(linmodel.darg.rho_c3(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.darg.rho_c3(ind.fch_d1:ind.fch_d2,:)))]);
% xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim([pav(ind.temp(1)), pav(ind.temp(2))]);
% setlabels(sprintf('[zoomed] LINEAR MODEL: PHASE DIFFERENCE OF \\rho_c TO f = %.1f MHz (IN DEGREE)', fch(ind.fch3)/1e6),...
%    'f [MHz]', 'P_{av} [W]', 'arg(\rho_c) [deg]');
% setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
% %     around settings.f4
% subplot(2,1,2); hold on;
% ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f4-settings.fd, 'nearest', 'extrap'));
% ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f4+settings.fd, 'nearest', 'extrap'));
% ind.temp = findzeros([1, all(isnan(arg.rho_c(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
% plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% surface(fch/1e6, pav, linmodel.darg.rho_c4', 'EdgeColor','none', 'FaceColor','interp');
% hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick,...
%    'CLim', [min(min(linmodel.darg.rho_c4(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.darg.rho_c4(ind.fch_d1:ind.fch_d2,:)))]);
% xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim([pav(ind.temp(1)), pav(ind.temp(2))]);
% setlabels(sprintf('[zoomed] LINEAR MODEL: PHASE DIFFERENCE OF \\rho_c TO f = %.1f MHz (IN DEGREE)', fch(ind.fch4)/1e6),...
%    'f [MHz]', 'P_{av} [W]', 'arg(\rho_c) [deg]');
% setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
% 
% % linear model: phase of difference rho_d vs f an Pav zoomed around some frequency values [1/2]
% figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear'];
% %     around settings.f1
% subplot(2,1,1); hold on;
% ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f1-settings.fd, 'nearest', 'extrap'));
% ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f1+settings.fd, 'nearest', 'extrap'));
% ind.temp = findzeros([1, all(isnan(arg.rho_d(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
% plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% surface(fch/1e6, pav, linmodel.darg.rho_d1', 'EdgeColor','none', 'FaceColor','interp');
% hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick);
% xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim([pav(ind.temp(1)), pav(ind.temp(2))]);
% setlabels(sprintf('[zoomed] LINEAR MODEL: PHASE DIFFERENCE OF \\rho_{\\Delta} TO f = %.1f MHz (IN DEGREE)', fch(ind.fch1)/1e6),...
%    'f [MHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [deg]');
% setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
% if ~isempty(settings.clim.darg_rho_d1)
%    set(gca, 'CLim', settings.clim.darg_rho_d1);
% else
%    set(gca, 'CLim', [min(min(linmodel.darg.rho_d1(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.darg.rho_d1(ind.fch_d1:ind.fch_d2,:)))]);
% end
% %     around settings.f2
% subplot(2,1,2); hold on;
% ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f2-settings.fd, 'nearest', 'extrap'));
% ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f2+settings.fd, 'nearest', 'extrap'));
% ind.temp = findzeros([1, all(isnan(arg.rho_d(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
% plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% surface(fch/1e6, pav, linmodel.darg.rho_d2', 'EdgeColor','none', 'FaceColor','interp');
% hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick);
% xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim([pav(ind.temp(1)), pav(ind.temp(2))]);
% setlabels(sprintf('[zoomed] LINEAR MODEL: PHASE DIFFERENCE OF \\rho_{\\Delta} TO f = %.1f MHz (IN DEGREE)', fch(ind.fch2)/1e6),...
%    'f [MHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [deg]');
% setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
% if ~isempty(settings.clim.darg_rho_d2)
%    set(gca, 'CLim', settings.clim.darg_rho_d2);
% else
%    set(gca, 'CLim', [min(min(linmodel.darg.rho_d2(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.darg.rho_d2(ind.fch_d1:ind.fch_d2,:)))]);
% end
% 
% % linear model: phase of difference rho_d vs f an Pav zoomed around some frequency values [2/2]
% figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear'];
% %     around settings.f3
% subplot(2,1,1); hold on;
% ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f3-settings.fd, 'nearest', 'extrap'));
% ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f3+settings.fd, 'nearest', 'extrap'));
% ind.temp = findzeros([1, all(isnan(arg.rho_d(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
% plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% surface(fch/1e6, pav, linmodel.darg.rho_d3', 'EdgeColor','none', 'FaceColor','interp');
% hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick);
% xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim([pav(ind.temp(1)), pav(ind.temp(2))]);
% setlabels(sprintf('[zoomed] LINEAR MODEL: PHASE DIFFERENCE OF \\rho_{\\Delta} TO f = %.1f MHz (IN DEGREE)', fch(ind.fch3)/1e6),...
%    'f [MHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [deg]');
% setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
% if ~isempty(settings.clim.darg_rho_d3)
%    set(gca, 'CLim', settings.clim.darg_rho_d3);
% else
%    set(gca, 'CLim', [min(min(linmodel.darg.rho_d3(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.darg.rho_d3(ind.fch_d1:ind.fch_d2,:)))]);
% end
% %     around settings.f4
% subplot(2,1,2); hold on;
% ind.fch_d1 = floor(interp1(fch, [1:1:settings.nfch], settings.f4-settings.fd, 'nearest', 'extrap'));
% ind.fch_d2 =  ceil(interp1(fch, [1:1:settings.nfch], settings.f4+settings.fd, 'nearest', 'extrap'));
% ind.temp = findzeros([1, all(isnan(arg.rho_d(ind.fch_d1:ind.fch_d2, :)), 1), 1] - 0.1)-1; % -0.1 => bias towards non-NaN
% plot3(fch/1e6, squeeze(pavopmin.all(:,  1)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% plot3(fch/1e6, squeeze(pavopmin.all(:,end)), 370*ones(settings.nfch,1), 'Color','k', 'LineWidth',1);
% surface(fch/1e6, pav, linmodel.darg.rho_d4', 'EdgeColor','none', 'FaceColor','interp');
% hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick);
% xlim([fch(ind.fch_d1), fch(ind.fch_d2)]/1e6); ylim([pav(ind.temp(1)), pav(ind.temp(2))]);
% setlabels(sprintf('[zoomed] LINEAR MODEL: PHASE DIFFERENCE OF \\rho_{\\Delta} TO f = %.1f MHz (IN DEGREE)', fch(ind.fch4)/1e6),...
%    'f [MHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [deg]');
% setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
% if ~isempty(settings.clim.darg_rho_d4)
%    set(gca, 'CLim', settings.clim.darg_rho_d4);
% else
%    set(gca, 'CLim', [min(min(linmodel.darg.rho_d4(ind.fch_d1:ind.fch_d2,:))), max(max(linmodel.darg.rho_d4(ind.fch_d1:ind.fch_d2,:)))]);
% end
