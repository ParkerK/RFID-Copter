% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% other - visualize setup of gate3d-simulations
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
% version = gate3d_vissetup()
%    Just returns the version number (string).
% gate3d_vissetup(readerpool, tagpool, loopsettings, surfaces, room_dim, visibility)
%    Creates a 3D plot of the current setup (with readers, tags, and surfaces). Figure window can be made
%    invisible by setting VISIBILITY='off'.
%
%
% ***** Interface definition *****
% function version = gate3d_vissetup(readerpool, tagpool, loopsettings, surfaces, room_dim, visibility)
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

function version = gate3d_vissetup(readerpool, tagpool, loopsettings, surfaces, room_dim, visibility)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% reader gain patterns and annotations
internalsettings.gain_diff       = 6; % degree gridlines
internalsettings.size_rdrpattern = 0.5; % m distance from center for normalized gain==1
internalsettings.txt_shift       = 0.6; % times maximum gain in x, y, z
internalsettings.tx_prefix       = ''; % prefix for annotation of transmitters
internalsettings.vtx_prefix      = ''; % prefix for annotation of virtual transmitters
 
% surfaces
internalsettings.surf_lightness  = [1, 0.9; 1000, 0.1]; % [n2 min, lightness; n2 max, lightness] surface lightness (the more reflective, the darker)
internalsettings.surf_gridres    = 0.5; % m surface gridline resolution (approximately)
internalsettings.surf_alpha      = 2/3; % alpha value (transparency 0..1) of surfaces

% misc
internalsettings.rdr_color =        get(0, 'DefaultAxesColorOrder'); % reader color order
internalsettings.rdr_color = [internalsettings.rdr_color; max(0,min(1,internalsettings.rdr_color*0.5))]; 
internalsettings.tag_color = [ones(10,1), 0.5*ones(10,1), zeros(10,1)]; % tag color order
internalsettings.annotate_txmax = 100; % maximum number of transmitters for annotation
internalsettings.fontsize       =   8; % font size for annotations 

% internalsettings.rdr_color = [0.1,0.1,0.3; 0.1,0.3,0.1; 0.3,0.3,0.1; 0.3,0.1,0.1];


% *******************************************************************************************************
% prepare

% reader antenna patterns
for i = 1 : readerpool.n
   if ~isempty(readerpool.ant{i})
      rdrant = load(readerpool.ant{i});
   else
      rdrant.az = 0 : internalsettings.gain_diff : 359;
      rdrant.el = 0 : internalsettings.gain_diff : 180;
      rdrant.gain = ones(length(rdrant.az), length(rdrant.el));
      rdrant.full3d = true;
   end
   if rdrant.full3d
      dp = round(internalsettings.gain_diff/mean(diff(rdrant.az)));
      [az, el] = meshgrid(rdrant.az(1:dp:end)+readerpool.ant_rot{i}(1), rdrant.el(1:dp:end)+readerpool.ant_rot{i}(2));
      [rdrpos_dir.x{i}, rdrpos_dir.y{i}, rdrpos_dir.z{i}] =...
         sph2cart(az*pi/180, pi/2-el*pi/180, rdrant.gain(1:dp:end, 1:dp:end)');
   else
      dp = round(internalsettings.gain_diff/mean(diff(rdrant.az)));
      rdrpos_dir.x{i} = (rdrant.el_gain(1:dp:end) .* sin((rdrant.el(1:dp:end)+readerpool.ant_rot{i}(2))*pi/180)) *...
         (rdrant.az_gain(1:dp:end) .* cos((rdrant.az(1:dp:end)+readerpool.ant_rot{i}(1))*pi/180))';
      rdrpos_dir.y{i} = (rdrant.el_gain(1:dp:end) .* sin((rdrant.el(1:dp:end)+readerpool.ant_rot{i}(2))*pi/180)) *...
         (rdrant.az_gain(1:dp:end) .* sin((rdrant.az(1:dp:end)+readerpool.ant_rot{i}(1))*pi/180))';
      rdrpos_dir.z{i} = (rdrant.el_gain(1:dp:end) .* cos((rdrant.el(1:dp:end)+readerpool.ant_rot{i}(2))*pi/180)) * (rdrant.az_gain(1:dp:end) )';
   end
   rdrpos_dir.c{i} = 2 * sqrt( rdrpos_dir.x{i}.^2 + rdrpos_dir.y{i}.^2 + rdrpos_dir.z{i}.^2 );
   rdrpos_dir.x{i} = internalsettings.size_rdrpattern * rdrpos_dir.x{i} + readerpool.pos{i}(1);
   rdrpos_dir.y{i} = internalsettings.size_rdrpattern * rdrpos_dir.y{i} + readerpool.pos{i}(2);
   rdrpos_dir.z{i} = internalsettings.size_rdrpattern * rdrpos_dir.z{i} + readerpool.pos{i}(3);
end

% tag positions
for i = 1 : tagpool.n
   tagpos.x{i} = loopsettings.pos_t{i}(1:loopsettings.n_sig, 1);
   tagpos.y{i} = loopsettings.pos_t{i}(1:loopsettings.n_sig, 2);
   tagpos.z{i} = loopsettings.pos_t{i}(1:loopsettings.n_sig, 3);
end

% surfaces (darker = more reflective)
surfpos.fieldnames    = fieldnames(surfaces);
surfpos.col_n2        = logspace(log10(internalsettings.surf_lightness(1,1)), log10(internalsettings.surf_lightness(2,1)), 128); 
surfpos.col_lightness = linspace(internalsettings.surf_lightness(1,2), internalsettings.surf_lightness(2,2), 128);
for i = 1 : length(surfpos.fieldnames)
   % lightness
   surfpos.lightness{i} = interp1(surfpos.col_n2, surfpos.col_lightness, surfaces.(surfpos.fieldnames{i}).n2, 'nearest', 'extrap');
   
   % position (roomsized if no bounds are given) [x; y; z] [min, max]
   %     copy bounds
   if isfield(surfaces.(surfpos.fieldnames{i}), 'bounds')
      bounds = surfaces.(surfpos.fieldnames{i}).bounds;
   else
      bounds = room_dim;
   end
   %     make sure the orientation is correct [x; y; z] [min, max]
   if size(bounds, 2) > size(bounds, 1)
      bounds = bounds';
   end
   %     add "0" as minimum if this is just a vector
   if size(bounds,2) == 1
      bounds = [zeros(3,1), bounds];
   end
   %     replace +/- Inf by room dimensions
   bounds(isinf(bounds)) = room_dim(isinf(bounds));
   %     center of the surface
   surfpos.center{i} = mean(bounds, 2);
   surfpos.center{i}(surfaces.(surfpos.fieldnames{i}).dim) = surfaces.(surfpos.fieldnames{i}).shift;

   switch surfaces.(surfpos.fieldnames{i}).dim
      % surface in yz
      case 1
         [surfpos.y{i}, surfpos.z{i}] = meshgrid(...
            linspace(bounds(2,1), bounds(2,2), max(2, ceil(diff(bounds(2,:))/internalsettings.surf_gridres))),...
            linspace(bounds(3,1), bounds(3,2), max(2, ceil(diff(bounds(3,:))/internalsettings.surf_gridres))));
         surfpos.x{i} = surfaces.(surfpos.fieldnames{i}).shift * ones(size(surfpos.y{i}));
         % surface in xz
      case 2
         [surfpos.x{i}, surfpos.z{i}] = meshgrid(...
            linspace(bounds(1,1), bounds(1,2), max(2, ceil(diff(bounds(1,:))/internalsettings.surf_gridres))),...
            linspace(bounds(3,1), bounds(3,2), max(2, ceil(diff(bounds(3,:))/internalsettings.surf_gridres))));
         surfpos.y{i} = surfaces.(surfpos.fieldnames{i}).shift * ones(size(surfpos.x{i}));
         % surface in xy
      case 3
         [surfpos.x{i}, surfpos.y{i}] = meshgrid(...
            linspace(bounds(1,1), bounds(1,2), max(2, ceil(diff(bounds(1,:))/internalsettings.surf_gridres))),...
            linspace(bounds(2,1), bounds(2,2), max(2, ceil(diff(bounds(2,:))/internalsettings.surf_gridres))));
         surfpos.z{i} = surfaces.(surfpos.fieldnames{i}).shift * ones(size(surfpos.x{i}));
      otherwise
         error('Unsupported surfaces.%s.dim = %i', surfpos.fieldnames{i}, surfaces.(surfpos.fieldnames{i}).dim);
   end
end



% *******************************************************************************************************
% create 3d plot

% figure
figure('Visible', visibility); hold on;

% readerpool
for i = 1 : readerpool.n
   % virtual transmitters: color similar to surfaces
   txt_shift = internalsettings.txt_shift;
   if readerpool.virt{i}
      mesh(rdrpos_dir.x{i}, rdrpos_dir.y{i}, rdrpos_dir.z{i}, 'EdgeColor',...
         max(0, min(1, internalsettings.rdr_color(readerpool.virt{i},:)+1-readerpool.virt_gf{i})));
      if readerpool.n < internalsettings.annotate_txmax
         text(readerpool.pos{i}(1) + txt_shift*internalsettings.size_rdrpattern, readerpool.pos{i}(2) - txt_shift*internalsettings.size_rdrpattern,...
            readerpool.pos{i}(3) + txt_shift*internalsettings.size_rdrpattern,...
            sprintf('%s%i', internalsettings.vtx_prefix, i), 'Fontsize',internalsettings.fontsize);
      end
   else
      mesh(rdrpos_dir.x{i}, rdrpos_dir.y{i}, rdrpos_dir.z{i}, rdrpos_dir.c{i}, 'EdgeColor', internalsettings.rdr_color(i,:));
      if readerpool.n < internalsettings.annotate_txmax
         text(readerpool.pos{i}(1) + txt_shift*internalsettings.size_rdrpattern, readerpool.pos{i}(2) - txt_shift*internalsettings.size_rdrpattern,...
            readerpool.pos{i}(3) + txt_shift*internalsettings.size_rdrpattern,...
            sprintf('%s%i', internalsettings.tx_prefix, i), 'Fontsize',internalsettings.fontsize);
      end
   end
end

% tag positions
for i = 1 : tagpool.n
   tag_color = internalsettings.tag_color(1+mod(i-1, size(internalsettings.tag_color,1)), :);
   plot3(tagpos.x{i}, tagpos.y{i}, tagpos.z{i}, '.', 'MarkerFaceColor',tag_color, 'MarkerEdgeColor',tag_color);
end

% surfaces (semi-transparent)
for i = 1 : length(surfpos.fieldnames)
   surface(surfpos.x{i}, surfpos.y{i}, surfpos.z{i}, 'EdgeColor','k',...
      'FaceColor', ones(3,1)*surfpos.lightness{i}, 'FaceAlpha', internalsettings.surf_alpha);
   if readerpool.n < internalsettings.annotate_txmax
      text(surfpos.center{i}(1), surfpos.center{i}(2),...
         surfpos.center{i}(3), surfpos.fieldnames{i}, 'Interpreter','none', 'Fontsize',internalsettings.fontsize);
   end
end

% annotation
setlabels('SIMULATION SETUP (VTX and surfaces: darker = higher gain)', 'x [m]', 'y [m]', 'z [m]');
hold off; grid on; view([40,30]); axis tight; axis equal;
lim.x = get(gca, 'Xlim'); lim.y = get(gca, 'Xlim'); lim.z = get(gca, 'Xlim');
xlim(zoom_out(get(gca, 'Xlim'))); ylim(zoom_out(get(gca, 'Ylim'))); zlim(zoom_out(get(gca, 'Zlim')));

% % plot circles
% center = zeros(1,3);
% for i = 1 : tagpool.n
%    center = center + [mean(tagpos.x{i}(:)), mean(tagpos.y{i}(:)), mean(tagpos.z{i}(:))] / tagpool.n;
% end
% hold on; plot3(center(1), center(2), center(3), 'k.');
% for radius = 5 : 5 : 100 % m
%    draw_circle_xz(center, radius, 'k:');
% end
% hold off;

end



% *******************************************************************************************************
% *******************************************************************************************************
% "zoom out" somewhat to better fit plot and text (works best for x,y,z in the "normal" directions)
function minmax = zoom_out(minmax)
   c =  sum(minmax)/2;
   d = diff(minmax)/2;
   minmax(1) = c - 1.05*d;
   minmax(2) = c + 1.15*d; % text is usually to the right side
end


% *******************************************************************************************************
% *******************************************************************************************************
% draw a circle in xz around a center
function draw_circle_xz(center, radius, linestyle)
   [x, z] = pol2cart(linspace(0, 2*pi, 250), radius);
   plot3(x+center(1), ones(size(x))+center(2), z+center(3), linestyle);
end
