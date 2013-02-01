% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% channel - location of readers/tags has changed => recalculate distances and setup channels accordingly
%
% References:
% [1] J.Karedal et al., "A Measurement-Based Statistical Model for Industrial Ultra-Wideband Channels,"
%     IEEE Tran Wireless Comm, vol 6, No 8, 2007, pp. 3028-3037
% [2] K. Witrisal, A New Method to Measure Parameters of Frequency-Selective Ration Channels Using Power
%     Measurements, IEEE Transactions on Communications, Vol. 59, No. 10, October 2001, pp. 1788-1800
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
% version = channel_newpos()
%    Just returns the version number (string).
% settings = channel_newpos(settings, pos_t, pos_r, ownsettings)
%    Modifies and returns the channel setup given in SETTINGS. Positions of transmitters and receivers
%    have to be provided in matrices/cell arrays POS_T and POS_R, respectively. Azimuth is set to 0 deg
%    and elevation is set to 90 deg by default (e.g., for 1-dim positions).
%    Smallscale settings are modified using linear interpolation and according to the vectors
%    given in OWNSETTINGS, for example:
%       ownsettings.v_dist = [ 0, 1, 2, 3]; ownsettings.v_k = [10, 9, 7, 0]; 
%       ownsettings.v_trms = [ 2, 4, 6, 8]; pos_t = {0, 1};  pos_r = {2, 3};  
%          => distance = [2, 1; 3, 2]
%          => Ricean K-factors = {7, 9; 0, 7};
%          => RMS delay spread = {6, 4; 8, 6};     
%    Reflective surfaces are added according to ownsettings.surfaces (if present). The reflected (NLOS) 
%    paths inherit all basic largescale and directivity properties (settings) of the direct (LOS) path.
%    
%
%
% ***** Interface definition *****
% function settings = channel_newpos(settings, pos_t, pos_r, ownsettings)
%    settings       channel settings to be modified (largescale, smallscale, directivity and surfaces)
%    pos_t          position of  N transmitters in meters; 
%                   e.g. 3-dim [x1, y1, z1; x2, y2, z2; ...; xN, yN, zN] or cell array of vectors 
%    pos_r          position of M receiver in meters; 
%                   e.g. 3-dim [x1, y1, z1; x2, y2, z2; ...; xN, yN, zN] or cell array of vectors 
%    ownsettings    struct containing settings for this function
%       .v_dist        distance reference vector in meters
%       .v_k           vector for Ricean K-factor "K = v_k(v_dist)" according to CHANNEL_SMALL
%       .v_trms        vector for RMS delay spread "trms = v_trms(v_dist)" according to CHANNEL_SMALL
%       .surfaces     (optional) struct containing substructs with definitions of one reflective surface each
%          .example       one such substruct, for example
%             .dim           dimension normal to reflective surface (e.g. 3 in [x,y,z] if xy reflects)
%             .shift         shift of this surface along the normal vector
%             .bounds        (optional) bounds for surface [min; max] for all dimensions
%                            ... dimension normal to the refl. surface will be ignored (zero width)
%
%    settings   modified largescale, smallscale, and reflection settings
%               (see CHANNEL_LARGE, CHANNEL_DIRECTIVITY, and CHANNEL_SURFACE)
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
% = make "VTX inherits trms of originating TX" optional (does not work for all scenarios)
% = simplify / split in subfunctions ... this function has become rather complex and confusing
% ? recalculate only modified positions (performance)
% ? vectorize calculation of POI
%
% *******************************************************************************************************


function settings = channel_newpos(settings, pos_t, pos_r, ownsettings)
version = 'beta 3.0';


% *******************************************************************************************************
% version system

% just return version number
if nargin == 0
   settings = version;
   return
end

% call version system (logs version numbers)
version_system();


% *******************************************************************************************************
% internal settings

% modes for surfaces (have to be identical to settings in CHANNEL_MAIN and CHANNEL_SURFACE)
internalsettings.surf_off  = 0; % neither reflection, nor transmission => surface inactive
internalsettings.surf_tran = 1; % "reflect" mode
internalsettings.surf_refl = 2; % "transmit" mode


% *******************************************************************************************************
% input parameter checks / prepare parameters

% backward compatibility: rename "reflections" to "surfaces"
if isfield(ownsettings, 'reflections')
   critwarn('Found obsolete field name "reflections" in ownsettings. Renamed to "surfaces". This check may be removed in the future.');
   ownsettings.surfaces = ownsettings.reflections;
   ownsettings = rmfield(ownsettings, 'reflections');   
end

% check contents of settings
%     prepare required data
expected.name = 'ownsettings';
expected.reqfields = {'v_dist', 'v_k', 'v_trms'};
%        reflective surfaces
if isfield(ownsettings, 'surfaces')
   expected.reqfields = [expected.reqfields, 'surfaces'];
   surf_names = fieldnames(ownsettings.surfaces);
   surf_cells = cell(0);
   for k = 1 : length(surf_names)
      surf_cells = [surf_cells, {surf_names{k}, {'dim', 'shift'}}];
   end
   expected.reqfields = [expected.reqfields, {'surfaces', surf_cells}];
   clear('surf_names', 'surf_cells');
end
%     check
errortext = contentcheck(ownsettings, expected);
%     output
if ~isempty(errortext)
   err('Incomplete settings\n%s', errortext);
end

% distance vector in ownsettings has to contain d=0 (otherwise saturation below won't work)
if min(ownsettings.v_dist) ~= 0  
   err('Distance vector in smallscale setup (call param. ownsettings) has to be positive at start at zero.');
end


% *******************************************************************************************************
% modify given setup

% cell -> array [x1,y1,z1; x2,y2,z2,; ...] if necessary to simplify calculations below
pos_t = pos_cell2mat(pos_t);
pos_r = pos_cell2mat(pos_r);

% calculate new distances and angles
[dist_tr, az_tr, el_tr] = get_distances(pos_t, pos_r);

% calculate distance for 1-point transmissions/reflections (if there is at least one surface)
if isfield(ownsettings, 'surfaces')
   % for all defined reflective surfaces
   surf_names = fieldnames(ownsettings.surfaces);
   % initialize
   poi = cellfun(@(x)  nan(size(pos_t, 1), size(pos_t, 2), size(pos_r, 1)), cell(length(surf_names),1), 'UniformOutput',false); % point of intersection
   moi =  cellfun(@(x) nan(size(pos_t, 1), size(pos_r, 1)), cell(length(surf_names),1), 'UniformOutput',false); % mode (type) of intersection
   for k = 1 : length(surf_names)
      % shortcuts
      dim_r = ownsettings.surfaces.(surf_names{k}).dim; % dimension of surface (index)
      shift = ownsettings.surfaces.(surf_names{k}).shift; % shift along the normal vector of this surface
      % mirrored positions
      %     transmitter
      pos_tm = pos_t;
      pos_tm(:, dim_r) = repmat(2*shift, size(pos_tm,1), 1) - pos_tm(:, dim_r);
      %     receiver
      pos_rm = pos_r;
      pos_rm(:, dim_r) = repmat(2*shift, size(pos_rm,1), 1) - pos_rm(:, dim_r);
      
      % get distances and angles
      [dist_tmr{k}, az_tmr{k}, el_tmr{k}] = get_distances(pos_tm, pos_r); %#ok<ASGLU,AGROW>
      [dist_trm{k}, az_trm{k}, el_trm{k}] = get_distances(pos_t, pos_rm); %#ok<ASGLU,AGROW>
      
      % get points of intersection
      for i = 1 : size(pos_t, 1)
         for j = 1 : size(pos_r, 1)
            % TX and RX in reflective plane
            if pos_t(i, dim_r) == shift && pos_r(j, dim_r) == shift
               warn('Transmitter %i and receiver %i in surface plane "%s". Switching off this surface.', i, j, surf_names{k});
               moi{k}(i, j)    = internalsettings.surf_off; % off
               poi{k}(i, :, j) = (pos_t(i, :) + pos_r(j, :)) / 2;
            % TX in reflective plane
            elseif pos_t(i, dim_r) == shift
               warn('Transmitter %i in surface plane "%s". Switching off this surface.', i, surf_names{k});
               moi{k}(i, j)    = internalsettings.surf_off; % off
               poi{k}(i, :, j) = pos_t(i, :);
            % RX in reflective plane
            elseif pos_r(j, dim_r) == shift
               warn('Receiver %i in surface plane "%s". Switching off this surface.', j, surf_names{k});
               moi{k}(i, j)    = internalsettings.surf_off; % off
               poi{k}(i, :, j) = pos_r(j, :);
            % otherwise: TX + unit vector * factor to reach intersection = poi
            else
               if sign(pos_t(i, dim_r) - shift) == sign(pos_r(j, dim_r) - shift) % RX and TX on the same side => point of reflection
                  moi{k}(i, j)    = internalsettings.surf_refl; % reflection
                  poi{k}(i, :, j) = pos_t(i, :) - pos_rm(j, :);
                  poi{k}(i, :, j) = poi{k}(i, :, j) / norm(poi{k}(i, :, j));
                  poi{k}(i, :, j) = pos_t(i, :) + poi{k}(i, :, j) * ( shift - pos_t(i, dim_r) ) /  poi{k}(i, dim_r, j);
               else % RX and TX on different sides of the surface => point of transmission
                  moi{k}(i, j)    = internalsettings.surf_tran; % transmission
                  poi{k}(i, :, j) = pos_t(i, :) - pos_r(j, :);
                  poi{k}(i, :, j) = poi{k}(i, :, j) / norm(poi{k}(i, :, j));
                  poi{k}(i, :, j) = pos_t(i, :) + poi{k}(i, :, j) * ( shift - pos_t(i, dim_r) ) /  poi{k}(i, dim_r, j);               
               end
            end
         end
  
      end
   end
   %     check: there should be no NaNs leftover in moi and poi
   if sum(sum( cell2mat(cellfun(@(x)      isnan(x), moi,  'UniformOutput',false)) )) +...
         sum(sum( cell2mat(cellfun(@(x) sum(isnan(x)), poi,  'UniformOutput',false)) )) ~= 0
      criterr('There are leftover NaNs in type/point of at least one intersection with a surface.');
   end
end

% warn if any distance is outside ownsettings.v_dist
if any(any(dist_tr > max(ownsettings.v_dist)))
   critwarn('Distances outside range given in ownsettings. Saturating trms, using linear extrap. for K.');
end

% for all combinations
for i = 1 : size(dist_tr, 1)
   for j = 1 : size(dist_tr, 2)
      
      % largescale and directivity model (if largescale model exists)
      if isfield(settings{i,j}, 'largescale')
         settings{i,j}.largescale.dist = dist_tr(i,j); % for all dimensions
         if ~isempty(az_tr)
            settings{i,j}.directivity.dir_tx(1) = az_tr(i,j) * 180/pi; % only for 2d,3d
         else
            settings{i,j}.directivity.dir_tx(1) = 0;
         end
         if ~isempty(el_tr)
            settings{i,j}.directivity.dir_tx(2) = el_tr(i,j) * 180/pi; % only for 3d
         else
            settings{i,j}.directivity.dir_tx(2) = 90; % default: in xy-plane
         end
      end
      
      % surfaces model, if exists (we also need largescale and directivity models for that)
      if isfield(ownsettings, 'surfaces') && isfield(settings{i,j}, 'largescale') && isfield(settings{i,j}, 'directivity')
         % process all surfaces for this connection TX -> RX
         surf_names = fieldnames(ownsettings.surfaces);
         for k = 1 : length(surf_names)
            
            % copy largescale and directivity settings
            settings{i,j}.surfaces.(surf_names{k}).largescale  = settings{i,j}.largescale;
            settings{i,j}.surfaces.(surf_names{k}).directivity = settings{i,j}.directivity;
            
            % copy surface settings
            settings{i,j}.surfaces.(surf_names{k}).dim    = ownsettings.surfaces.(surf_names{k}).dim;
            
            % modify settings
            if moi{k}(i,j) == internalsettings.surf_refl; % reflection (TX->RXm and RX->TXm)
               [settings{i,j}.surfaces.(surf_names{k}).largescale, settings{i,j}.surfaces.(surf_names{k}).directivity] = ...
                  surf_modsettings(dist_tmr{k}(i,j), az_tmr{k}(i,j), az_trm{k}(i,j), el_tmr{k}(i,j), el_trm{k}(i,j),...
                  settings{i,j}.surfaces.(surf_names{k}).largescale, settings{i,j}.surfaces.(surf_names{k}).directivity);
            else % transmission (TX->RX and RX->TX) or offline 
               [settings{i,j}.surfaces.(surf_names{k}).largescale, settings{i,j}.surfaces.(surf_names{k}).directivity] = ...
                  surf_modsettings(dist_tr(i,j), az_tr(i,j), az_tr(i,j), el_tr(i,j), el_tr(i,j),...
                  settings{i,j}.surfaces.(surf_names{k}).largescale, settings{i,j}.surfaces.(surf_names{k}).directivity);
            end
            
            % calculate angle of incidence on surface
            switch settings{i,j}.surfaces.(surf_names{k}).dim
               case 1 % "x => yz-plane"
                  alpha = settings{i,j}.surfaces.(surf_names{k}).directivity.dir_tx(1);
               case 2 % "y => xz-plane"
                  alpha = settings{i,j}.surfaces.(surf_names{k}).directivity.dir_tx(1) - 90;
               case 3 % "z => xy-plane"
                  alpha = settings{i,j}.surfaces.(surf_names{k}).directivity.dir_tx(2);
               otherwise
                  err('Unsupported normal vector dimension "%i" for surface.', settings{i,j}.surfaces.(surf_names{k}).dim);
            end
            settings{i,j}.surfaces.(surf_names{k}).aoi = abs(round(alpha/180)*180 - alpha);
            
            % calculate distance to the shortest edge and set mode to transmission or reflection
            settings{i,j}.surfaces.(surf_names{k}).poi2e = poi2edge(ownsettings.surfaces.(surf_names{k}), poi{k}(i,:,j));
            settings{i,j}.surfaces.(surf_names{k}).mode = moi{k}(i,j); % mode of intersection: refl. or trans. (if surface was active)
            
         end
      end
      
      % smallscale model
      settings{i,j}.smallscale.k = interp1(ownsettings.v_dist, ownsettings.v_k,...
         dist_tr(i,j), 'linear', 'extrap');
      %     VTX/VRX: inherit RMS delay spread from originating TX/RX (will create roughly constant decay for low K)
      if isfield(settings{i,j}, 'virtual') && (settings{i,j}.virtual.vtx || settings{i,j}.virtual.vrx)
         if settings{i,j}.virtual.vtx
            settings{i,j}.smallscale.trms = interp1(ownsettings.v_dist, ownsettings.v_trms,...
               dist_tr(settings{i,j}.virtual.otx,j), 'linear', ownsettings.v_trms(end));
         end
         if settings{i,j}.virtual.vrx
            settings{i,j}.smallscale.trms = interp1(ownsettings.v_dist, ownsettings.v_trms,...
               dist_tr(i,settings{i,j}.virtual.orx), 'linear', ownsettings.v_trms(end));
         end
      else
         settings{i,j}.smallscale.trms = interp1(ownsettings.v_dist, ownsettings.v_trms,...
            dist_tr(i,j), 'linear', ownsettings.v_trms(end));
      end
   end
end

end



% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% convert cell array of positions to matrix of positions ([x1,y1,z1; x2,y2,z2,; ...])
%    ... not pretty, but robust
function m = pos_cell2mat(c)
if iscell(c)
   m = zeros(length(c), length(c{1}));
   for i = 1 : length(c)
      m(i,:) = c{i}(:)';
   end
else
   m = c;
end
end



% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% modify channel settings for surface: directivity
%     dist: distance TX*-RX*; az_* azimuth angle [rad]; el_* elevation angle [rad]
%     *_tmr: mirrored TX, original RX; *_trm: original TX, mirrored RX, 
%     FOR REFLECTION:   surf_modsettings(dist_tr,  az_tr,  az_tr,  el_tr,  el_tr,  largescale, directivity)
%     FOR TRANSMISSION: surf_modsettings(dist_trm, az_tmr, az_trm, el_tmr, el_trm, largescale, directivity)
function [largescale, directivity] = surf_modsettings(dist, az_tmr, az_trm, el_tmr, el_trm, largescale, directivity)
largescale.dist = dist;
if ~isempty(az_tmr) % 2d,3d
   directivity.dir_tx(1) = az_trm * 180/pi; % TX looks at the mirrored image of RX
   directivity.dir_rx(1) = az_tmr * 180/pi; % ... and vice-versa
else
   directivity.dir_tx(1) = 0;
end
if ~isempty(el_tmr) % 3d
   directivity.dir_tx(2) = el_trm * 180/pi; % TX looks at the mirrored image of RX
   directivity.dir_rx(2) = el_tmr * 180/pi; % ... and vice-versa
else
   directivity.dir_tx(2) = 90; % default: in xy-plane
end
end



% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% check point of intersection: returns shortest distance of the intersection to the closest edge
%   surface: ownsettings struct of this surface
%   poi:     point of intersection
function poi2e = poi2edge(surface, poi)
% surface has no bounds: there is definitely a POI
if ~isfield(surface, 'bounds')
   poi2e = -Inf; % really far inside the surface
   return
end
% otherwise: check if the POI is within all dimensions of the surface
if isfield(surface, 'bounds')
   %     all the relevant dimensions
   dim = 1 : numel(poi);
   dim(surface.dim) = [];
   %     get all relevant distances to surface borders (negative: inside surface)
   dist = [poi(dim)-surface.bounds(2, dim), surface.bounds(1, dim)-poi(dim)];
   %     get distance to the closest edge
   if all( dist <= 0 ) % if inside the surface
      poi2e = max(dist);
   else
      poi2e = min(dist(dist>0));   
   end
end
end



% *******************************************************************************************************
% *******************************************************************************************************
% DEBUG

% PLOT OF THE ENTIRE SETUP (INCLUDTING VTX AND VRX)
%       if isfield(ownsettings.surfaces.(surf_names{k}), 'bounds');
%          b = ownsettings.surfaces.(surf_names{k}).bounds;
%       else
%          b = [-10,10; -10,10; -10,10];
%       end
%       switch dim_r
%          case 1
%             [y, z] = meshgrid(linspace(b(1,2),b(2,2),11), linspace(b(1,3),b(2,3),11));
%             x = zeros(size(y)) + shift;
%          case 2
%             [x, z] = meshgrid(linspace(b(1,1),b(2,1),11), linspace(b(1,3),b(2,3),11));
%             y = zeros(size(x)) + shift;
%          case 3
%             [x, y] = meshgrid(linspace(b(1,1),b(2,1),11), linspace(b(1,2),b(2,2),11));
%             z = zeros(size(x)) + shift;
%       end
%       figure(1); clf(1); hold on;
%       %       figure; hold on;
%       surface(x, y, z, 'FaceColor', 'none'); % reflective plane
%       for i = 1 : size(pos_t, 1)
%          for j = 1 : size(pos_r, 1)
%             plot3(pos_t(i,1), pos_t(i,2), pos_t(i,3), 'bsquare'); % transmitters
%             plot3(pos_r(j,1), pos_r(j,2), pos_r(j,3), 'bo'); % receivers
%             plot3([pos_t(i,1), pos_r(j,1)], [pos_t(i,2), pos_r(j,2)], [pos_t(i,3), pos_r(j,3)], 'g-'); % TX -> RX
%             plot3([pos_t(i,1), pos_rm(j,1)], [pos_t(i,2), pos_rm(j,2)], [pos_t(i,3), pos_rm(j,3)], 'b-'); % TX -> RXm
%             plot3([pos_tm(i,1), pos_r(j,1)], [pos_tm(i,2), pos_r(j,2)], [pos_tm(i,3), pos_r(j,3)], 'c-'); % TXm -> RX
%             plot3(poi{k}(i, 1, j), poi{k}(i, 2, j), poi{k}(i, 3, j), 'ro');
%          end
%       end
%       hold off; xlabel('x'); ylabel('y'); zlabel('z'); title(sprintf('%s', upper(surf_names{k})));
%       legend('surface', 'TX', 'RX', 'TX -> RX', 'TX -> RXm', 'TXm -> RX', 'intersection');
%       view([25,50]);
%       pause



% *******************************************************************************************************
% *******************************************************************************************************
% OLD AND RUSTY

