% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% demo analysis script, IEEE RFID 2012
%
% *******************************************************************************************************

% *******************************************************************************************************
% initialization
clear; close all; clc; pause(0.01);

% paths, etc
%     path to globalinit.m
%pathstring = '~/paris-osf';
pathstring = '/Users/johnpohollaren/Documents/Reynolds_Lab/RFID-Copter/PARIS/source/';

%     add this path and call globalinit script
restoredefaultpath; addpath(pathstring); clear('dummy', 'path');
globalsettings_copy = globalinit('silent');

% "switch to here"
cd(fileparts(mfilename('fullpath')));



% *******************************************************************************************************
% load data and reshape to 3D representation

% load data
load /Users/johnpohollaren/Documents/Reynolds_Lab/RFID-Copter/PARIS/source/simulator/_data/mfcw_gate3d/john

% receiver positions
for i = 1 : settings.tagpool.n
   rx_pos{i}.x  = settings.loop.pos_t{i}(:, 1);
   rx_pos{i}.y  = settings.loop.pos_t{i}(:, 2);
   rx_pos{i}.z  = settings.loop.pos_t{i}(:, 3);
   rx_pos{i}.az = settings.tagpool.ant_rot{i}(1) * ones(settings.loop.n, 1);
   rx_pos{i}.el = settings.tagpool.ant_rot{i}(2) * ones(settings.loop.n, 1);
end

% incident power level at the tag (narrowband)
ind_rdr = find(cellfun(@(x) x==0, settings.readerpool.virt)); % skip virtual transmitters
tag_pwr = nan(length(ind_rdr), settings.tagpool.n, settings.loop.n);
for r = 1 : length(ind_rdr)
   for t = 1 : settings.tagpool.n
      for l = 1 : settings.loop.n
         tag_pwr_dbm(r,t,l) = 10*log10(settings.readerpool.ptx{r}) + 30 + ...
            20*log10(abs(sum( results.channel_rt{r,t}.cir_g{l} .*...
            exp(complex(0, 2*pi*settings.readerpool.fc{r} * results.channel_rt{r,t}.cir_d{l})) )));
      end
   end
end

% reshape to 2D arrays supported by surface()
for t = 1 : settings.tagpool.n
   nx{t} = length(unique(rx_pos{t}.x));
   ny{t} = length(unique(rx_pos{t}.y));
   nz{t} = length(unique(rx_pos{t}.z));
%    if t == 1
   if nx{t} == 1
      x{t} = reshape(rx_pos{t}.x, nz{t}, ny{t});
      y{t} = reshape(rx_pos{t}.y, nz{t}, ny{t});
      z{t} = reshape(rx_pos{t}.z, nz{t}, ny{t});
      for r = 1 : length(ind_rdr)
         p{r,t} = reshape(squeeze(tag_pwr_dbm(r,t,:))', nz{t}, ny{t});
         trms{r,t} = reshape(results.channel_rt{r,t}.trms', nz{t}, ny{t}) * 1e9;
      end
   elseif ny{t} == 1
      x{t} = reshape(rx_pos{t}.x, nz{t}, nx{t});
      y{t} = reshape(rx_pos{t}.y, nz{t}, nx{t});
      z{t} = reshape(rx_pos{t}.z, nz{t}, nx{t});
      for r = 1 : length(ind_rdr)
         p{r,t} = reshape(squeeze(tag_pwr_dbm(r,t,:))', nz{t}, nx{t});
         trms{r,t} = reshape(results.channel_rt{r,t}.trms', nz{t}, nx{t}) * 1e9;
      end
%    elseif t == 2
   elseif nz{t} == 1
      x{t} = reshape(rx_pos{t}.x, ny{t}, nx{t});
      y{t} = reshape(rx_pos{t}.y, ny{t}, nx{t});
      z{t} = reshape(rx_pos{t}.z, ny{t}, nx{t});
      for r = 1 : length(ind_rdr)
         p{r,t} = reshape(squeeze(tag_pwr_dbm(r,t,:))', ny{t}, nx{t});
         trms{r,t} = reshape(results.channel_rt{r,t}.trms', ny{t}, nx{t}) * 1e9;
      end
   end
   pavg{t} = 10.^(p{1,t}/10) / length(ind_rdr);
   for r = 2 : length(ind_rdr)
      pavg{t} = pavg{t} + 10.^(p{r,t}/10) / length(ind_rdr);
   end
   pavg{t} = 10*log10(pavg{t});
end

% interpolate
usf = 4;     
for t = 1 : settings.tagpool.n
   % basis
   [u1, u2] = meshgrid(1:size(x{t}, 1), 1:size(x{t}, 2));
   [v1, v2] = meshgrid(1:1/usf:size(x{t}, 1), 1:1/usf:size(x{t}, 2));
   % positions
   xi{t} = interp2(u1,u2, x{t}.', v1, v2, 'spline').';
   yi{t} = interp2(u1,u2, y{t}.', v1, v2, 'spline').';
   zi{t} = interp2(u1,u2, z{t}.', v1, v2, 'spline').';
   % data
   for r = 1 : length(ind_rdr)
      pi{r,t} = interp2(u1,u2, p{r,t}.', v1, v2, 'spline').';
      trmsi{r,t} = interp2(u1,u2, trms{r,t}.', v1, v2, 'spline').';
   end
   pavgi{t} = interp2(u1,u2, pavg{t}.', v1, v2, 'spline').';
end


% return

% *******************************************************************************************************
% output / plots

% number of positions, simulation time
fprintf('\nSimulation time: %.1f minutes for %d positions (%d reader-tag connections; %.0f ms per connection)\n', results.misc.simtime/60, ...
   sum(cellfun(@length, settings.loop.pos_t)), sum(cellfun(@length, settings.loop.pos_t)) * length(ind_rdr),...
   1e3 * results.misc.simtime / ( sum(cellfun(@length, settings.loop.pos_t)) * length(ind_rdr) ) );

close all
% 3D plot with setup
settings.loop.n = 0; % remove tags
settings.loop.n_sig = 0;
gate3d_vissetup(settings.readerpool, settings.tagpool, settings.loop, settings.channel_global.surfaces, settings.room_dim, 'on');
hold on;
for t = 1 : settings.tagpool.n
   surface(xi{t}, yi{t}, zi{t}, pavgi{t}, 'EdgeColor','none', 'FaceColor', 'interp');
end
hold off; setcolorbar(); set(gca, 'clim', [-15,10]);
setlabels('AVG TAG RX POWER [dBM] FOR ALL READERS','x [m]',' y [m]',' z[m]');
view([-50,20]);  saveas(gcf, 'mfcw_gate3d__live-rfid12_3d_tagpwr', 'fig');
zlim([0 5])

