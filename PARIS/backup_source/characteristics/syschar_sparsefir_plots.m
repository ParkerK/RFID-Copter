% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% generate plots for: system - sparse FIR implementation
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


% ******************************************************************************************************
% initialization
clear; close all; clc;
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
globalsettings_copy = globalinit('exceptions');

% "switch to here"
cd(fileparts(mfilename('fullpath')));


% ******************************************************************************************************
% settings

plotsettings.hostname     = 'fredcolon'; % for which host?
plotsettings.folder       = 'plots'; % subfolder for plots
plotsettings.lenf_slices  = [10,20,50,100,200,500,1000]; % place slices at these filter lengths (exec. time plots) 

% color limits, to make plots between hosts comparable (leave empty for automatic)
plotsettings.clim_log     = [-6, 3]; % log10 seconds execution time


% ******************************************************************************************************
% load data and create plots

% load workspace
load(fullfile('workspaces', sprintf('syschar_sparsefir_%s_workspace', plotsettings.hostname)));

% tipping point plot
figure; 
surface(v_len_i, v_len_f, results.sxing', 'FaceColor','interp');
axis('tight'); set(gca, 'xscale','log', 'yscale','log', 'CLim',[0,1]);
setlabels('SPARSITY AT PERFORMANCE CROSSING FILTER/LOOP', 'input signal length [samples]', 'filter length [samples]',...
   'sparsity @ crossing [%] (below: filter faster)'); setcolorbar();
savefigure(gcf, fullfile(plotsettings.folder, sprintf('syschar_sparsefir_%s_sxing', plotsettings.hostname)) , 'png')

% execution time plots
%     reformat to 3d structure
for i = 1 : settings.n_i
   for f = 1 : settings.n_f
      results.t_filter(i,f,:) = 10.^interp1(results.loop{i,f}.sparsity, log10(results.loop{i,f}.t_filter), s_int, 'linear', 'extrap');
      results.t_loop(i,f,:)   = 10.^interp1(results.loop{i,f}.sparsity, log10(results.loop{i,f}.t_loop),   s_int, 'linear', 'extrap');
   end
end
%     get logarithmic color limits
if isempty(plotsettings.clim_log)
   clim_log = [min(min(min(log10([results.t_filter, results.t_loop])))), ...
               max(max(max(log10([results.t_filter, results.t_loop]))))];
else
   clim_log = plotsettings.clim_log;
end
%     plot (filter)
figure;
hx = slice(len_f, len_i, s_int*100, log10(results.t_filter), plotsettings.lenf_slices, [], []);
set(hx, 'FaceColor','interp','EdgeColor','none'); set(gca, 'xscale','log','yscale','log', 'clim',clim_log);
axis tight; view(20, 60); setcolorbar();
setlabels('LOGARITHMIC EXECUTION TIME FOR FILTER IMPLEMENTATION (extrap!)', 'filter length [samples]',...
   'input signal length [samples]', 'sparsity [%]');
savefigure(gcf, fullfile(plotsettings.folder, sprintf('syschar_sparsefir_%s_t-filter', plotsettings.hostname)) , 'png');
%     plot (loop)
figure;
hx = slice(len_f, len_i, s_int*100, log10(results.t_loop), plotsettings.lenf_slices, [], []);
set(hx, 'FaceColor','interp','EdgeColor','none'); set(gca, 'xscale','log','yscale','log', 'clim',clim_log);
axis tight; view(20, 60); setcolorbar();
setlabels('LOGARITHMIC EXECUTION TIME FOR LOOP IMPLEMENTATION (extrap!)', 'filter length [samples]',...
   'input signal length [samples]', 'sparsity [%]');
savefigure(gcf, fullfile(plotsettings.folder, sprintf('syschar_sparsefir_%s_t-loop', plotsettings.hostname)) , 'png');

% cleanup
% close all;
