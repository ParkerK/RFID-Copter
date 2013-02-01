% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate plots for: tag - modulation "light" characteristic (reflection coefficient)
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
settings.version = 'beta 3.0';

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


% *******************************************************************************************************
% settings

%     general / file settings
settings.plotfolder   = 'plots'; % subfolder for plots (there will be additional subfolders for the plots)
settings.reset_folder = true; % delete folder before saving new plots?
settings.charid       = '_c-p2a_d'; % start with '_' if not empty c-p2a_d
%     list of assembly combinations to plot (combination cat(1),rat(1) ... cat(2),rat(2) ... plotted)
settings.char_ind.rat = [-98, -95, -90, -80, -67, -50, -33, 0, 50, 100, 200, 500, 1000, 2000, 5000, 10000];
settings.char_ind.cat = 450*ones(size(settings.char_ind.rat));
%     list of detuning combinations to plot (combination cat(1),rat(1) ... cat(2),rat(2) ... plotted)
settings.char_ind.enr = -0.9 : 0.1 : 1; % enhancement of resonance
settings.char_ind.fsr =  0 : 10e6 : 200e6; % Hz frequency shift of resonance

% MFCW distance error estimates
settings.mfcw.df12 = 1e6; % Hz frequency differences between 2FCW carrier pair

% pav min/max and tick
settings.pav_min = 1e-5; % W
settings.pav_max =    3; % W
settings.pav_tick  = [1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1e0, 1e1]; % W

% other
settings.rfid_range = [864, 954]*1e6; % Hz RFID frequency range
%     blank out phases close to zero magnitude (will jump)
settings.arg_magmin = 1e-2; % min. magnitude before phase is replaced by NaN

% order of unwrapping for phase, leave empty to deactivate unwrapping
settings.unwrap_rho = [2,1,3]; % reflection coefficient (1: f, 2: mod, 3: power); default: [2,1,3]
settings.unwrap_lin =   [1,2]; % linear model (1: f, 2: power); default: [1,2]

% clim for some plots (can be empty: [] for auto)
settings.clim.rho_c_mag  = [  0,    1]; % rho_c vs f and pav (linear model)
settings.clim.rho_d_mag  = [  0, 0.40]; % rho_d vs f and pav (linear model)
settings.clim.rho_c_arg  = [-pi,   2*pi]; % rad rho_c vs f and pav (linear model)
settings.clim.rho_d_arg  = [-pi, 1.5*pi]; % rad rho_d vs f and pav (linear model)   
settings.clim.mfcw_derr1 = [-2.0, 2.0]; % m 2FCW distance error ("full span" plot)
settings.clim.mfcw_derr2 = [-0.4, 0.4]; % m 2FCW distance error (zoomed plot)
settings.clim.darg_rho   = [ -10,  10]; % deg zoomed phase difference: rho_d(mod) - rho_d(unmod)
% 
% settings.clim.rho_c_mag  = []; % rho_c vs f and pav (linear model)
% settings.clim.rho_d_mag  = []; % rho_d vs f and pav (linear model)
% settings.clim.rho_c_arg  = []; % rad rho_c vs f and pav (linear model)
% settings.clim.rho_d_arg  = []; % rad rho_d vs f and pav (linear model)   
% settings.clim.mfcw_derr1 = []; % m 2FCW distance error ("full span" plot)
% settings.clim.mfcw_derr2 = []; % m 2FCW distance error (zoomed plot)
% settings.clim.darg_rho   = []; % deg zoomed phase difference: rho_d(mod) - rho_d(unmod)

% figure options
settings.figopt.colormap        = jet(128); % used colormap
settings.figopt.visible         =     'on'; % figure visibility {'on', 'off'}
% settings.figopt.papersize       =  [14, 7]; % cm (for papers)
settings.figopt.type            =     'png'; % type of figures to create ('png' or 'eps'; leave empty to create both; see SAVEFIGURE)
settings.figopt.nonconfidential =     false; % remove confidential information (e.g. for presentations)?
settings.figopt.resolution      =       150; % dpi image resolution
% settings.figopt.addspace_xlabel = 0;
% settings.figopt.addspace_ylabel = 0;


% SETTINGS FOR THESIS
settings.charid       = '_c-p2a_d'; % start with '_' if not empty c-p2a_d
settings.char_ind.fsr =    [0e6, 75e6]; % Hz frequency shift of resonance
settings.char_ind.enr =    [0.0,  0.5]; % enhancement of resonance
settings.char_ind.rat =    [0,0]; % assembly resistance shift
settings.char_ind.cat = 250e-15*[1,1]; % assemby capacity
settings.pav_min      = 1e-5; % W
settings.pav_max      =    1; % W
settings.figopt.papersize       = [14.5, 10]; % cm size [width, height] 
settings.figopt.visible         = 'off'; % figure visibility {'on', 'off'} 
settings.figopt.type            = 'eps'; % type of figures to create ('png' or 'eps'; leave empty to create both; see SAVEFIGURE)
settings.figopt.nonconfidential = false; % remove confidential information (e.g. for presentations)?
settings.figopt.resolution      =   250; % dpi image resolution
settings.figopt.addspace_ylabel =  0.00;
settings.figopt.addspace_colorb =  0.03;



% *******************************************************************************************************
% load light characteristic and process data

% output own name/version first
disp(sprintf('*******************************************************************************************************'));
disp(sprintf('* This is %s, version %s', mfilename, settings.version));
disp(sprintf('*******************************************************************************************************'));

% load characteristic (this might take some time)
disp(sprintf('\n Loading data (this might take a while)...'));
load(sprintf('tagchar_modulator%s-light', settings.charid));
%     split file: ignore
if exist('splitfile', 'var')
   tagchar_mod = tagchar_mod_chunk;
end

% loop for set characteristics
for loop_ind = 1 : min([length(settings.char_ind.cat), length(settings.char_ind.rat), length(settings.char_ind.enr), length(settings.char_ind.fsr)])
   
   % map assembly and detuning to states
%    settings.state_a = find( round(ads_mapping.assembly.v_cat*1e15) == settings.char_ind.cat(loop_ind) &...
%       ads_mapping.assembly.v_rats == settings.char_ind.rat(loop_ind) );
%    settings.state_d = find( abs(ads_mapping.detuning.v_res_en - settings.char_ind.enr(loop_ind)) < 1e-10 &...
%       ads_mapping.detuning.v_fshift == settings.char_ind.fsr(loop_ind) );
   [settings.state_a, settings.state_d] = adstate('settings->state', ads_mapping,...
      settings.char_ind.cat(loop_ind), settings.char_ind.rat(loop_ind),...
      settings.char_ind.enr(loop_ind), settings.char_ind.fsr(loop_ind));
   
   %     check
   if isempty(settings.state_a)
      error('Assembly impedance not part of the state vectors.');
   end
   if isempty(settings.state_d)
      error('Detuning state for not of the state vectors.');
   end
   
   % map assembly/detuning state to index in characteristic
   char_ind = find( sweepsettings.a_full == settings.state_a & sweepsettings.d_full == settings.state_d );
   %     check
   if isempty(char_ind)
      error('Could not find this assembly/detunign state in the characteristic.');
   end
   
   % extract data and output what is done
   fch = tagchar_mod{char_ind}.fch;
   pav = tagchar_mod{char_ind}.pav;
   pic = tagchar_mod{char_ind}.pic;
   disp(sprintf('\n\n*******************************************************************************************************'));
   disp(sprintf('* WORKING ON: %s.mat [LIGHT] (created by %s)', tagchar_mod{char_ind}.matfilename, tagchar_mod{char_ind}.createdby));
   
   
   % ************************************************************
   % prepare data: linear model (rho_c + rho_{delta})
   disp(sprintf('   Preparing data for linear model...'));
   
   % decompose to linear model
   rho_c = ( squeeze(tagchar_mod{char_ind}.rho_pav(:, 1, :)) + squeeze(tagchar_mod{char_ind}.rho_pav(:, end, :)) ) / 2;
   rho_d =   squeeze(tagchar_mod{char_ind}.rho_pav(:, 1, :)) - rho_c;
   
   % magnitue
   mag.rho     = abs(tagchar_mod{char_ind}.rho_pav);
   mag.rho_c   = abs(rho_c);
   mag.rho_d   = abs(rho_d);
   
   % phase
   arg.rho     = angle(tagchar_mod{char_ind}.rho_pav);
   arg.rho_c   = angle(rho_c);
   arg.rho_d   = angle(rho_d);
   %     unwrap
   for i = settings.unwrap_rho
      arg.rho    (end:-1:1,:,end:-1:1) = unwrap(arg.rho    (end:-1:1,:,end:-1:1), [], i);
   end
   for i = settings.unwrap_lin
      arg.rho_c  (end:-1:1,  end:-1:1) = unwrap(arg.rho_c  (end:-1:1,  end:-1:1), [], i);
      arg.rho_d  (end:-1:1,  end:-1:1) = unwrap(arg.rho_d  (end:-1:1,  end:-1:1), [], i);
   end
   
   % blank out phases of linear models when magnitude is close to zero
   arg.rho_c( mag.rho_c < settings.arg_magmin ) = NaN;
   arg.rho_d( mag.rho_d < settings.arg_magmin ) = NaN;
   
   % differentiate
   fch_drho = (fch(1:end-1) + fch(2:end)) / 2;
   dfch = diff(fch);
   mag.drho_c = nan(tagchar_mod{char_ind}.settings.nfch-1, tagchar_mod{char_ind}.settings.np);
   arg.drho_c = nan(tagchar_mod{char_ind}.settings.nfch-1, tagchar_mod{char_ind}.settings.np);
   mag.drho_d = nan(tagchar_mod{char_ind}.settings.nfch-1, tagchar_mod{char_ind}.settings.np);
   arg.drho_d = nan(tagchar_mod{char_ind}.settings.nfch-1, tagchar_mod{char_ind}.settings.np);
   for j = 1 : tagchar_mod{char_ind}.settings.np
      mag.drho_c(:, j) = diff(mag.rho_c(:, j)) ./ dfch;
      arg.drho_c(:, j) = diff(arg.rho_c(:, j)) ./ dfch;
      mag.drho_d(:, j) = diff(mag.rho_d(:, j)) ./ dfch;
      arg.drho_d(:, j) = diff(arg.rho_d(:, j)) ./ dfch;
   end
   
   
   % ************************************************************
   % prepare data: minimum available power levels
   disp(sprintf('   Calculating minimum available power necessary for tag operation...'));
   
   pavopmin.all      = zeros(tagchar_mod{char_ind}.settings.nfch, tagchar_mod{char_ind}.settings.nm);
   pavopmin.square50 = zeros(tagchar_mod{char_ind}.settings.nfch,           1);
   for i = 1 : tagchar_mod{char_ind}.settings.nfch
      % Pav,min vs. frequency and state of modulation
      for j = 1 : tagchar_mod{char_ind}.settings.nm
         pavopmin.all(i,j) = interp1(pic, squeeze(tagchar_mod{char_ind}.m_pav(i, j, :)), tagchar_mod{char_ind}.picopmin(i), 'linear', 'extrap');
      end
   end
   
   
   % ************************************************************
   % MFCW distance errors based on linear model
   % d_err = (phi_d0-phi_d1) * c/(2*w_d);
   disp(sprintf('   Estimating MFCW distance errors based on linear model...'));
   
   % 2FCW distance error; fixed frequency difference vs. center frequency (f0, f0+df)
   %     prepare cell array of matrices
   for k = 1 : length(settings.mfcw.df12)
      linmodel.mfcw.derr_f0{k} = nan(tagchar_mod{char_ind}.settings.nfch, tagchar_mod{char_ind}.settings.np);
   end
   %     calculate
   warning('off', 'MATLAB:interp1:NaNinY');
   warning('off', 'MATLAB:chckxy:IgnoreNaN');
   for i = 1 : tagchar_mod{char_ind}.settings.nfch
      for j = 1 : tagchar_mod{char_ind}.settings.np
         for k = 1 : length(settings.mfcw.df12)
            linmodel.mfcw.derr_f0{k}(i, j) = 3e8./(4*pi*settings.mfcw.df12(k)) * ...
               ( arg.rho_d(i, j) - interp1(fch, arg.rho_d(:, j), fch(i)+settings.mfcw.df12(k), 'linear','extrap') );
         end
      end
   end
   warning('on', 'MATLAB:interp1:NaNinY');
   warning('on', 'MATLAB:chckxy:IgnoreNaN');
   
   


%    close all
% 
% 
%    % linear model: angle of difference rho_d vs f an Pav (for papers)
%    figure('Visible', settings.figopt.visible, 'Colormap',settings.figopt.colormap); hold on; plotnames = {'linear_diff-phase'};
% %    plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
% %    plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
%    surface(fch/1e9, pav, 180/pi*arg.rho_d', 'EdgeColor','none', 'FaceColor','interp');
%    grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim([min(fch),max(fch)]/1e9); ylim([settings.pav_min, settings.pav_max]);
%    setlabels('', 'f [GHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [deg]');
%    setlegend({'P_{min} mod, unmod'}, 'SouthEast'); setcolorbar();
%    
%       % save and close figures
%    %   create new directory, meanwhile deactivate warning "Removed '/blah' from the MATLAB path..."
%    warning('off', 'MATLAB:RMDIR:RemovedFromPath');
%    [status, message, messageid] = rmdir(fullfile(settings.plotfolder, sprintf('tagchar_mod%s_light', tagchar_mod{char_ind}.settings.suffix)), 's'); %#ok<NASGU>
%    [status, message, messageid] = mkdir(fullfile(settings.plotfolder, sprintf('tagchar_mod%s_light', tagchar_mod{char_ind}.settings.suffix)));
%    warning('on', 'MATLAB:RMDIR:RemovedFromPath');
%    %   save figures
%    for i = 1 : gcf
%       if isempty(plotnames{i})
%          savefigure(i, fullfile(settings.plotfolder, sprintf('tagchar_mod%s_light', tagchar_mod{char_ind}.settings.suffix),...
%             sprintf('tagchar_mod%s_light_%02i', tagchar_mod{char_ind}.settings.suffix, i)), settings.figopt.type, settings.figopt);
%       else
%          savefigure(i, fullfile(settings.plotfolder, sprintf('tagchar_mod%s_light', tagchar_mod{char_ind}.settings.suffix),...
%             sprintf('tagchar_mod%s_light_%02i_%s', tagchar_mod{char_ind}.settings.suffix, i, plotnames{i})),...
%             settings.figopt.type, settings.figopt);
%       end
% %       close(i);
%    end
%    
   
  
   % ************************************************************
   % plots
   disp(sprintf('   Creating plots...')); plotnames = {}; 
   warning('off', 'MATLAB:legend:UnsupportedFaceColor');
      
   % rho vs f an Pav for Zmod=min and Zmod=max (non-extrapolated)
   %     magnitude
   figure('Visible', settings.figopt.visible, 'Colormap',settings.figopt.colormap); plotnames = [plotnames, 'rho_nonext-mag'];
   subplot(2,1,1); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, squeeze(mag.rho(:, end, :))', 'EdgeColor','interp', 'FaceColor','interp');
   hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]); setcolorbar; %view([80,45]);
   set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'cLim',xyzlimits(mag.rho(:, [1,end], :)));
   setlabels('MAGNITUDE OF UNMODULATED REFLECTION COEFFICIENT','f [GHz]', 'P_{av} [W]', '|\rho|');
   setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   subplot(2,1,2); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, squeeze(mag.rho(:,   1, :))', 'EdgeColor','interp', 'FaceColor','interp');
   hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]); setcolorbar; %view([80,45]);
   set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'cLim',xyzlimits(mag.rho(:, [1,end], :)));
   setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   setlabels('MAGNITUDE OF MODULATED REFLECTION COEFFICIENT','f [GHz]', 'P_{av} [W]', '|\rho|');
   %     phase
   figure('Visible', settings.figopt.visible, 'Colormap',settings.figopt.colormap); plotnames = [plotnames, 'rho_nonext-phase'];
   subplot(2,1,1); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 2.1*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 2.1*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, squeeze(arg.rho(:, end, :))', 'EdgeColor','interp', 'FaceColor','interp');
   hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]); setcolorbar; %view([80,45]);
   set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'cLim',xyzlimits(arg.rho(:, [1,end], :)));
   setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   setlabels('PHASE OF UNMODULATED REFLECTION COEFFICIENT IN RAD','f [GHz]', 'P_{av} [W]', 'arg(\rho) [rad]');
   subplot(2,1,2); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 2.1*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 2.1*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, squeeze(arg.rho(:,   1, :))', 'EdgeColor','interp', 'FaceColor','interp');
   hold off; grid on; xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]); setcolorbar; %view([80,45]);
   set(gca, 'YScale','log', 'YTick',settings.pav_tick, 'cLim',xyzlimits(arg.rho(:, [1,end], :)));
   setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   setlabels('PHASE OF MODULATED REFLECTION COEFFICIENT IN RAD','f [GHz]', 'P_{av} [W]', 'arg(\rho) [rad]');
   
   % linear model:center rho_c vs f an Pav with fixed color limits for comparison
   figure('Visible', settings.figopt.visible, 'Colormap',settings.figopt.colormap); plotnames = [plotnames, 'linear_center-comp'];
   %     magnitude
   subplot(2,1,1); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, mag.rho_c', 'EdgeColor','none', 'FaceColor','interp');
   grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
   setlabels('LINEAR MODEL: CENTER VALUE \rho_c (MAGNITUDE) FOR COMPARISON', 'f [GHz]', 'P_{av} [W]', '|\rho_c|');
   setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   if ~isempty(settings.clim.rho_c_mag); set(gca, 'CLim', settings.clim.rho_c_mag); end;
   %     phase
   subplot(2,1,2); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 2.1*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 2.1*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, arg.rho_c', 'EdgeColor','none', 'FaceColor','interp');
   grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
   setlabels('LINEAR MODEL: CENTER VALUE \rho_c (PHASE IN RAD) FOR COMPARISON', 'f [GHz]', 'P_{av} [W]', 'arg(\rho_c) [rad]');
   setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   if ~isempty(settings.clim.rho_c_arg); set(gca, 'CLim', settings.clim.rho_c_arg); end;
   
   % linear model: difference rho_d vs f an Pav with fixed color limits for comparison
   figure('Visible', settings.figopt.visible, 'Colormap',settings.figopt.colormap); plotnames = [plotnames, 'linear_diff-comp'];
   %     magnitude
   subplot(2,1,1); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 1.1*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 1.1*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, mag.rho_d', 'EdgeColor','none', 'FaceColor','interp');
   grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
   setlabels('LINEAR MODEL: DIFFERENCE \rho_{\Delta} (MAGNITUDE) FOR COMPARISON', 'f [GHz]', 'P_{av} [W]', '|\rho_{\Delta}|');
   setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   if ~isempty(settings.clim.rho_d_mag); set(gca, 'CLim', settings.clim.rho_d_mag); end;
   %     phase
   subplot(2,1,2); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 2.1*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 2.1*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, arg.rho_d', 'EdgeColor','none', 'FaceColor','interp');
   grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
   setlabels('LINEAR MODEL: DIFFERENCE \rho_{\Delta} (PHASE IN RAD) FOR COMPARISON', 'f [GHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [rad]');
   setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   if ~isempty(settings.clim.rho_d_arg); set(gca, 'CLim', settings.clim.rho_d_arg); end;
   
% % %    % linear model: angle of difference rho_d vs f an Pav (for papers)
% % %    figure('Visible', settings.figopt.visible, 'Colormap',settings.figopt.colormap); hold on; plotnames = [plotnames, 'linear_diff-phase'];
% % %    plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
% % %    plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*pi*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
% % %    surface(fch/1e9, pav, 180/pi*arg.rho_d', 'EdgeColor','none', 'FaceColor','interp');
% % %    grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
% % %    setlabels('', 'f [GHz]', 'P_{av} [W]', 'arg(\rho_{\Delta}) [deg]');
% % %    setlegend({'P_{min} mod, unmod'}, 'NorthEast'); setcolorbar();
   
   disp('MANUAL MODIFICATIONS')
   % linear model: separate plot for drho_d (phase)
   figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_diff-phase']; hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(settings.rfid_range(1)*[1,1]/1e9, [settings.pav_min,settings.pav_max], 370*[1,1], 'k--',  'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(settings.rfid_range(2)*[1,1]/1e9, [settings.pav_min,settings.pav_max], 370*[1,1], 'k--',  'LineWidth',1);
   surface(fch/1e9, pav, arg.rho_d'*180/pi, 'EdgeColor','none', 'FaceColor','interp');
   hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
   setlabels('', 'f_c [GHz]', 'P_{av} [W]'); setlegend({'P_{av,min} (mod, unmod)', 'RFID freq. range'}, 'NorthEast'); setcolorbar();
   
   disp('MANUAL MODIFICATIONS')
   % linear model: separate plot for drho_d/df (phase)
   figure('Visible', settings.figopt.visible); plotnames = [plotnames, 'linear_grad-phase']; hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 370*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(settings.rfid_range(1)*[1,1]/1e9, [settings.pav_min,settings.pav_max], 370*[1,1], 'k--',  'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 370*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(settings.rfid_range(2)*[1,1]/1e9, [settings.pav_min,settings.pav_max], 370*[1,1], 'k--',  'LineWidth',1);
   surface(fch_drho/1e9, pav, arg.drho_d'*1e6*180/pi, 'EdgeColor','none',  'FaceColor','interp');
   hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
   setlabels('', 'f_c [GHz]', 'P_{av} [W]'); setlegend({'P_{av,min} (mod, unmod)', 'RFID freq. range'}, 'NorthEast'); setcolorbar();
   
   % MFCW: distance error for carrier at fch, one secondary carrier at fch+df12 (2FCW)
   figure('Visible', settings.figopt.visible, 'Colormap',settings.figopt.colormap); plotnames = [plotnames, 'mfcw_hist-derr-f0'];
   %     overview 
   subplot(2,1,1); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 110*max(max(linmodel.mfcw.derr_f0{1}))*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 110*max(max(linmodel.mfcw.derr_f0{1}))*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, 100*linmodel.mfcw.derr_f0{1}', 'EdgeColor','none', 'FaceColor','interp');
   hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
   setlabels(sprintf('LINEAR MODEL: DISTANCE ERROR [cm] FOR 2FCW, FREQ. DIFFERENCE f_1-f_0 = %.1f MHz', settings.mfcw.df12(1)/1e6),...
      'f_0 [GHz] (main carrier)', 'P_{av} [W]', 'd_{err} [cm]');
   setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
   if ~isempty(settings.clim.mfcw_derr1); set(gca, 'CLim', 100*settings.clim.mfcw_derr1); end;
   %     zoomed color limits
   subplot(2,1,2); hold on;
   plot3(fch/1e9, squeeze(pavopmin.all(:,  1)), 110*max(max(linmodel.mfcw.derr_f0{1}))*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   plot3(fch/1e9, squeeze(pavopmin.all(:,end)), 110*max(max(linmodel.mfcw.derr_f0{1}))*ones(tagchar_mod{char_ind}.settings.nfch,1), 'Color','k', 'LineWidth',1);
   surface(fch/1e9, pav, 100*linmodel.mfcw.derr_f0{1}', 'EdgeColor','none', 'FaceColor','interp');
   hold off; grid on; set(gca, 'YScale','log', 'YTick',settings.pav_tick); xlim(xyzlimits(fch)/1e9); ylim([settings.pav_min, settings.pav_max]);
   setlabels('[ZOOMED COLOR LIMITS]', 'f_0 [GHz] (main carrier)', 'P_{av} [W]', 'd_{err} [cm]');
   setlegend({'P_{min} mod, unmod'}, 'SouthWest'); setcolorbar();
   if ~isempty(settings.clim.mfcw_derr2); set(gca, 'CLim', 100*settings.clim.mfcw_derr2); end;
   
   
   % ************************************************************
   % save plots to file
   disp(sprintf('   Saving plots...'));
   warning('on', 'MATLAB:legend:UnsupportedFaceColor');
   
   % check if all plots have been named
   if gcf ~= length(plotnames)
      error('Mismatch between number of figures and length of plotnames array.');
   end
   
   % save and close figures
   %   create new directory, meanwhile deactivate warning "Removed '/blah' from the MATLAB path..."
   warning('off', 'MATLAB:RMDIR:RemovedFromPath');
   [status, message, messageid] = rmdir(fullfile(settings.plotfolder, sprintf('tagchar_mod%s_light', tagchar_mod{char_ind}.settings.suffix)), 's'); %#ok<NASGU>
   [status, message, messageid] = mkdir(fullfile(settings.plotfolder, sprintf('tagchar_mod%s_light', tagchar_mod{char_ind}.settings.suffix)));
   warning('on', 'MATLAB:RMDIR:RemovedFromPath');
   %   save figures
   for i = 1 : gcf
      if isempty(plotnames{i})
         savefigure(i, fullfile(settings.plotfolder, sprintf('tagchar_mod%s_light', tagchar_mod{char_ind}.settings.suffix),...
            sprintf('tagchar_mod%s_light_%02i', tagchar_mod{char_ind}.settings.suffix, i)), settings.figopt.type, settings.figopt);
      else
         savefigure(i, fullfile(settings.plotfolder, sprintf('tagchar_mod%s_light', tagchar_mod{char_ind}.settings.suffix),...
            sprintf('tagchar_mod%s_light_%02i_%s', tagchar_mod{char_ind}.settings.suffix, i, plotnames{i})),...
            settings.figopt.type, settings.figopt);
      end
      close(i);
   end
   
end

