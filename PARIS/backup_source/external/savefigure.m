% *******************************************************************************************************
% Position Aware RFID Systems: The PARIS Simulation Framework
% ***********************************************************
% resizes a figure, filling as much space as possible, and saves the figure to different formats
%
% Important:
%   = Do not maximize figure window before calling this function. Maximized windows cannot be resized.
%   = No sanity checks are done; inconsistent settings may lead to strange behavior.
%   = Axes and legends/colorbars have to be created subsequently ("zipper principle").   
%   = If colorbars/legends are used in subplots, each subplot has to have a colorbar/legend.
%   = For colorbars, only location "EastOutside" is supported.
%   = The function can be deactivated via a switch: optional.off=true
%
% Comments:
%   - This function may have different behavior on different operating systems, because saved/printed 
%     figures are NOT 100% identical to displayed figures. The smaller the spacings, the larger the
%     chance that labels are outside the saved image.
%   - Reported positions may be completely wrong after resize operations (more likely on Linux?). There
%     are some workarounds and checks in place to deal with this issue.
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
% savefigure(fig, filename, filetype, optional)
%    Resizes figure FIG (maximize used space) and saves the result as graphic FILETYPE (optional; 
%    default: eps plus png) in file FILENAME. Set FILETYPE to '' in order to be able to use OPTIONAL and
%    generate the default filetypes.
%
%
% ***** Interface definition *****
% function savefigure(fig, filename, filetype, optional)
%    fig             figure handle
%    filename        filename for graphics
%    filetype        (optional) string {'png', 'pdf', 'jpg', 'eps', 'epsbw', '+fig'=png+eps+fig} default: png+eps
%    optional        (optional) struct with custom settings; also see internalsettings in this function
%       .papersize          figure size in cm (default: [25, 15.4] cm)
%       .resolution         image resolution in dpi (default: 300 dpi)
%       .axisgap            empty space around each plot axis+colobar [left, bottom, right, top]
%                           (normalized units; default: [0.000, 0.010, 0.010, 0.002])
%       .axisgap3d          empty space around each plot axis+colorbar [left, bottom, right, top]
%                           (normalized units; default: [0.000, 0.040, 0.010, 0.002])
%       .colorbarwidth      size of colorbars (normalized units; default: 0.02)
%       .maxcolorbarwidth   maximum size of colorbars in cm (default: 3 cm)
%       .colorbargap        empty space between axes and colorbars (normalized units; default: 0.005)
%       .legendaddwidth     additional width for legends (normalized units; default: 0.001)
%       .legendgap          empty space between legend and axis (normalized units; default: 0.01)
%       .font               replace all fonts in figure with this font, if present (default: no replacement)
%       .waittimefactor     additional wait time to make sure the figure is ready after resizing it 
%                           (times tic...toc for resize operation; default: 2)
%       .waittime           maximum additional wait time according to .waittimefactor in sec (default: 0.2)
%       .nonconfidential    if true, all tick labels, colorbars and titles will be removed (default: false)
%       .off                return without doing anything (naturally false per default)
%
%
% ***** Changelog *****
% REVISION   DATE         USER        DESCRIPTION (! bugfixes, + addons, - removals, ~ otherwise)
% beta 1.0   2010-03-30   arnitz      ~ initial release
% beta 2.0   2010-09-01   arnitz      ~ testing release (unstable)
% beta 3.0   2012-05-07   arnitz      ~ partial bugfix release
% beta 3.1   2012-05-22   arnitz      ~ complete revision of the entire function
%                                     + optional override for all settings via OPTIONAL, additional settings
%                                     + support for external png->eps conversion (way better results)
% beta 3.2   2012-06-15   arnitz      ! colorbar/legend handles are made lowercase (again) before search;
%                                       something seems to change the case, and Tag property is case sen.
% beta 3.3   2012-06-26   arnitz      ~ eps conversion also for eps+png and +fig type
%
%
% ***** Todo *****
% ! internalsettings.legendaddwidth has influence on horz legend position. Why?!
% = detect legend/colorbar alignment to be able to support arbitrary combinations of legends/colorbars
% - shift label text closer to axis in 3D plots (possible?)
% ? x and y label outside outerposition for for small papersizes
%
% *******************************************************************************************************

function savefigure(fig, filename, filetype, optional)

% *******************************************************************************************************
% internal settings (Matlab figure -> graphic is unfortunatly not WYSIWYG)

% size / alignment
%     figure
internalsettings.papersize  = [25, 15.4]; % paper size cm size [width, height] 
internalsettings.resolution =        300; % image resolution [dpi]
%     axes / subplots (times figure width); will be scaled by aspect ratio
internalsettings.axisgap   = [0.000, 0.010, 0.010, 0.002]; % space around each axis [left, bottom, right, top]
internalsettings.axisgap3d = [0.000, 0.040, 0.010, 0.002]; % space around each axis [left, bottom, right, top] for 3D plots
%     colorbar
internalsettings.colorbarwidth  =  0.02; % size of colorbars (times figure width)
internalsettings.maxcolbarwidth =     3; % cm maximum size of colorbar (useful for huge figures)
internalsettings.colorbargap    = 0.005; % gap between colorbar and axes border (times figure width)
%     legends
internalsettings.legendgap      = 0.01; % gap between legend and axes border (times figure width)
internalsettings.legendaddwidth = 0.001; % additional width for legends (times figure width)

% alignment / checks
%     alignment
internalsettings.snapdf_legend = 0.67; % snap distance factor for legends (in height and width dimension)
%     checks
internalsettings.sizetol    = 0.05; % size tolerance for resize operation (relative error)
internalsettings.tilepostol = 0.075; % tolerance for position in subplot tiling detection

% misc
%     figure (re-)size timing problem workarounds: wait for figure to be complete
internalsettings.waittimefactor =   2; % wait times the time it took for the set(fig, ...resize) command
internalsettings.maxwaittime    = 0.2; % s maximum wait time after each resize operation
%     output settings
internalsettings.jpegqual = 92; % jpeg quality (if saved as jpg); 0..100
%     system/conversion commands
internalsettings.cmd_png2eps     = 'convert -compress zip %s.png eps2:%s.eps'; % command for 
internalsettings.cmd_png2eps_chk = 'convert --version > /dev/null'; % command to check availability of this tool
internalsettings.cmd_png2eps_rm  = 'rm -f %s.png'; % command for removing the temporary png file
%     misc
internalsettings.nonconfidential = false; % remove "confidential" information (tick labels, title, colorbars) if true
internalsettings.off = false; % return without doing anything (can be used to globally switch of saving of figures)



% *******************************************************************************************************
% input parameter checks / prepare parameters

if nargin < 2
   error('Not enough input arguments');
elseif nargin == 2
   filetype = '';
elseif nargin == 4
   internalsettings = mergestructs(internalsettings, optional); % externally controlled settings
end

% add (limited) backward compatibility
internalsettings = backward_compatibility(internalsettings);

% switched off?
if internalsettings.off
   return;
end

% complete/modify internalsettings
%     aspect ratio (make relative x and y distances identical in created figure)
internalsettings.aspectratio = internalsettings.papersize(1) / internalsettings.papersize(2);
%     scale gaps with aspect ratio to create identical gaps in final plot
internalsettings.axisgap([1,3]) = internalsettings.axisgap([1,3]) / sqrt(internalsettings.aspectratio);
internalsettings.axisgap([2,4]) = internalsettings.axisgap([2,4]) * sqrt(internalsettings.aspectratio);
internalsettings.axisgap3d([1,3]) = internalsettings.axisgap3d([1,3]) / sqrt(internalsettings.aspectratio);
internalsettings.axisgap3d([2,4]) = internalsettings.axisgap3d([2,4]) * sqrt(internalsettings.aspectratio);
internalsettings.colorbargap = internalsettings.colorbargap * [1/sqrt(internalsettings.aspectratio), sqrt(internalsettings.aspectratio)];
internalsettings.legendgap   = internalsettings.legendgap   * [1/sqrt(internalsettings.aspectratio), sqrt(internalsettings.aspectratio)];
%     saturate colorbar size
internalsettings.colorbarwidth = min(internalsettings.colorbarwidth, internalsettings.maxcolbarwidth/internalsettings.papersize(1));

% font change requested?
if isfield(internalsettings, 'font');
   setfonts(gcf, internalsettings.font);
end

% would it be a good idea to use an external conversion program?
internalsettings.conv_png = false;
if ~any(strcmpi(filetype, {'png', 'pdf', 'jpg', 'epsbw'}))  &&...
      ~strcmpi(get(fig, 'renderer'), 'painters') % ... only 'painters' can produce vector outputs
   % check if png->eps conversion tool is available; if so, use this tool instead of Matlab's implementation
   if ~system(internalsettings.cmd_png2eps_chk) % return value 0: everyting OK
      internalsettings.conv_png = true;
      if strcmpi(filetype, 'eps')
         filetype = 'png';
      end
   else
      warning('savefigure:conversion', 'Direct EPS output of this figure will likely produce a large file. An attempt to use an external conversion tool failed.');
   end
end


% *******************************************************************************************************
% collect handles

% make sure we can see all handles
old.root.shh = get(0, 'ShowHiddenHandles');
set(0, 'ShowHiddenHandles', 'on');

% get all Axes objects in fig (this will include colobars and legends)
handles.all = findobj(fig, 'Type', 'Axes');

% get legends and colorbars (have to be tagged beforehand)
%     make all tags lowercase, just to make sure (Tag property is case sensitive; Type property is not) 
for i = 1 : length(handles.all)
   set(handles.all(i), 'Tag', lower(get(handles.all(i), 'Tag')));
end
%     get handles
handles.legends   = findobj(fig, 'Tag', 'legend');
handles.colorbars = findobj(fig, 'Tag', 'colorbar');

% make legends and colorbars invisible for axes search (otherwise legends and colorbars would register as axes)
set(0, 'ShowHiddenHandles', 'off');
setall(handles.legends,   'HandleVisibility', 'off');
setall(handles.colorbars, 'HandleVisibility', 'off');

% search for axes
handles.axes = findobj(fig, 'Type', 'axes');

% check number of ... compared to axes
%     ... colorbars
if ~isempty(handles.colorbars) && length(handles.axes) > length(handles.colorbars)
   error('More axes than colorbars (and at least one colorbar): not supported.')
end
%     ... legends
if ~isempty(handles.legends) && length(handles.axes) > length(handles.legends)
   error('More axes than legends (and at least one legend): not supported.')
end

% restore root visibility property to its previous state
set(0, 'ShowHiddenHandles', old.root.shh);



% *******************************************************************************************************
% remove "confidential" information if requested
if internalsettings.nonconfidential
   % remove colorbars
   for i = 1 : length(handles.colorbars)
      delete(handles.colorbars(i));
   end
   handles.colorbars = [];
   % remove titles and tick labels
   for i = 1 : length(handles.axes)
      delete(get(handles.axes(i), 'Title'));
      set(handles.axes(i), 'XTickLabel','', 'YTickLabel','', 'ZTickLabel','');
   end
   % redraw
   pause(internalsettings.waittime);
end



% *******************************************************************************************************
% resize figure

% make sure all necessary properties of all objects are set properly
for i = 1 : length(handles.all)
   set(handles.all(i), 'ActivePositionProperty','OuterPosition', 'Units','normalized', 'HandleVisibility','on');
end

% resize operation
ticID = tic;
set(fig,...
   'ActivePositionProperty', 'OuterPosition',...
   'Units', 'centimeters',...
   'Position', [0,0,internalsettings.papersize],... % +1 to avoid screen border
   'PaperUnits', 'centimeters',...
   'PaperPositionMode', 'manual',...
   'PaperSize', internalsettings.papersize,...
   'PaperPosition', [0,0,internalsettings.papersize]);
drawnow;
internalsettings.waittime = toc(ticID) * internalsettings.waittimefactor;

% make sure figure is ready ... otherwise the reported sizes (position, outerposition, ...) will be wrong
pause(min(internalsettings.maxwaittime, internalsettings.waittime));

% check if size is okay
pos = get(fig, 'Position');
if abs(pos(3)-internalsettings.papersize(1))/internalsettings.papersize(1) > internalsettings.sizetol ||...
      abs(pos(4)-internalsettings.papersize(2))/internalsettings.papersize(2) > internalsettings.sizetol
   error(['Unable to set figure size. Possible problems: figure window maximized, screen resolution too small,',...
      ' timing problems. Try increasing waittimefactor and maxwaittime if you suspect timing issues.']);
end



% *******************************************************************************************************
% collect layout information

% OuterPosition [left, bottom, width, height]
outerposition.axes      = getall(handles.axes,      'OuterPosition', 4);
outerposition.colorbars = getall(handles.colorbars, 'OuterPosition', 4);
outerposition.legends   = getall(handles.legends,   'OuterPosition', 4);

% Position [left, bottom, width, height]
position.axes      = getall(handles.axes,      'Position', 4);
position.colorbars = getall(handles.colorbars, 'Position', 4);
position.legends   = getall(handles.legends,   'Position', 4);

% TightInset [left, bottom, right, top]
tightinset.axes      = getall(handles.axes,      'TightInset', 4);
tightinset.colorbars = getall(handles.colorbars, 'TightInset', 4);
tightinset.legends   = getall(handles.legends,   'TightInset', 4);

% View [azimuth, elevation}
view.axes = getall(handles.axes, 'View', 2);

% discern number of subplots
tiles.horz = length(unique_tol( outerposition.axes(:,1), internalsettings.tilepostol ));
tiles.vert = length(unique_tol( outerposition.axes(:,2), internalsettings.tilepostol ));
%     check
if length(handles.axes) > tiles.horz * tiles.vert
   error('Detection of subplot alignment (horz, vert) failed.');
end

% check if colorbar position is supported
if ~isempty(outerposition.colorbars) && ...
      any(outerposition.colorbars(:, 3) > outerposition.colorbars(:, 4)) % width > hight?
   error('Positioning of colorbar not supported ... select "EastOutside".');
end



% *******************************************************************************************************
% calculate new outerpositions and positions

% prepare arrays
%     outerpositions
newouterposition.axes = zeros(length(handles.axes), 4);
newouterposition.colorbars = zeros(length(handles.colorbars), 4);
newouterposition.legends = zeros(length(handles.legends), 4);
%     positions
newposition.axes = zeros(length(handles.axes), 4);
newposition.colorbars = zeros(length(handles.colorbars), 4);
newposition.legends = zeros(length(handles.legends), 4);

% outerposition of axes
for i = 1: length(handles.axes)
   if prod(view.axes(i,:))~=0 % 3d plot: different spacing (need more space for labels)
      % from left and lower corner
      newouterposition.axes(i, 1) = round( outerposition.axes(i, 1) * tiles.horz ) / tiles.horz + internalsettings.axisgap3d(1);
      newouterposition.axes(i, 2) = round( outerposition.axes(i, 2) * tiles.vert ) / tiles.vert + internalsettings.axisgap3d(2);
      % width and height
      newouterposition.axes(i, 3) = round( outerposition.axes(i, 3) * tiles.horz ) / tiles.horz - internalsettings.axisgap3d(1) - internalsettings.axisgap3d(3);
      newouterposition.axes(i, 4) = round( outerposition.axes(i, 4) * tiles.vert ) / tiles.vert - internalsettings.axisgap3d(2) - internalsettings.axisgap3d(4);
   else
      % from left and lower corner
      newouterposition.axes(i, 1) = round( outerposition.axes(i, 1) * tiles.horz ) / tiles.horz + internalsettings.axisgap(1);
      newouterposition.axes(i, 2) = round( outerposition.axes(i, 2) * tiles.vert ) / tiles.vert + internalsettings.axisgap(2);
      % width and height
      newouterposition.axes(i, 3) = round( outerposition.axes(i, 3) * tiles.horz ) / tiles.horz - internalsettings.axisgap(1) - internalsettings.axisgap(3);
      newouterposition.axes(i, 4) = round( outerposition.axes(i, 4) * tiles.vert ) / tiles.vert - internalsettings.axisgap(2) - internalsettings.axisgap(4);
   end
end

% horizontal alignment of colorbars and axes (assumes that axis and colorbar handles are aligned) 
for i = 1: length(handles.colorbars)
   % width of colorbar (width-tightinset = textwidth)
   newouterposition.colorbars(i, 3) = internalsettings.colorbarwidth + outerposition.colorbars(i, 3) - tightinset.colorbars(i, 1) - tightinset.colorbars(i, 3) ;
   % shift to right border of subplot
   newouterposition.colorbars(i, 1) = newouterposition.axes(i, 1) + newouterposition.axes(i, 3) - newouterposition.colorbars(i, 3);
   % correct width of axis
   newouterposition.axes(i, 3) = newouterposition.axes(i, 3) - newouterposition.colorbars(i, 3) - internalsettings.colorbargap(1);
end

% position of axes
for i = 1: length(handles.axes)
   newposition.axes(i, 1) = newouterposition.axes(i, 1) + tightinset.axes(i, 1);
   newposition.axes(i, 2) = newouterposition.axes(i, 2) + tightinset.axes(i, 2);
   newposition.axes(i, 3) = newouterposition.axes(i, 3) - tightinset.axes(i, 1) - tightinset.axes(i, 3);
   newposition.axes(i, 4) = newouterposition.axes(i, 4) - tightinset.axes(i, 2) - tightinset.axes(i, 4);
end

% final layout of colorbars
for i = 1: length(handles.colorbars)
   % outerposition: vertical arrangement identical to position axis
   newouterposition.colorbars(i, 2) = newposition.axes(i, 2);
   newouterposition.colorbars(i, 4) = newposition.axes(i, 4);
   % position
   newposition.colorbars(i, 1) = newouterposition.colorbars(i, 1);
   newposition.colorbars(i, 2) = newouterposition.colorbars(i, 2);
   newposition.colorbars(i, 3) = newouterposition.colorbars(i, 3) - tightinset.colorbars(i, 3); % labels
   newposition.colorbars(i, 4) = newouterposition.colorbars(i, 4);
end

% position and layout of legends (assumes that axis and legend handles are aligned)
for i = 1: length(handles.legends)
   % width, height: use current
   newouterposition.legends(i, 3) = outerposition.legends(i, 3) + internalsettings.legendaddwidth;
   newouterposition.legends(i, 4) = outerposition.legends(i, 4);
   
   % position: shift to closest alignment position ({north/south, east/west})
   %     determine distance of legend borders to axis borders x=[left, right], y=[bottom, top]
   distx = [outerposition.legends(i, 1) - position.axes(i, 1), position.axes(i, 1) + position.axes(i, 3) - outerposition.legends(i, 1) - outerposition.legends(i, 3)];
   disty = [outerposition.legends(i, 2) - position.axes(i, 2), position.axes(i, 2) + position.axes(i, 4) - outerposition.legends(i, 2) - outerposition.legends(i, 4)];
   %     final position: x
   if distx(1) <= distx(2) * internalsettings.snapdf_legend % much closer to left side: place west
      newouterposition.legends(i, 1) = newposition.axes(i, 1) + internalsettings.legendgap(1);
   elseif distx(2) <= distx(1) * internalsettings.snapdf_legend % much closer to right side: place east
      newouterposition.legends(i, 1) = newposition.axes(i, 1) + newposition.axes(i, 3) - newouterposition.legends(i, 3) - internalsettings.legendgap(1);
   else % distance to both borders approximately equal: place in center
      newouterposition.legends(i, 1) = newposition.axes(i, 1) + newposition.axes(i, 3)/2 - newouterposition.legends(i, 3)/2;
   end
   %     final position: y
   if disty(1) <= disty(2) * internalsettings.snapdf_legend % much closer to bottom: place south
      newouterposition.legends(i, 2) = newposition.axes(i, 2) + internalsettings.legendgap(2);
   elseif disty(2) <= disty(1) * internalsettings.snapdf_legend % much closer to top: place north
      newouterposition.legends(i, 2) = newposition.axes(i, 2) + newposition.axes(i, 4) - newouterposition.legends(i, 4) - internalsettings.legendgap(2);
   else % distance to both borders approximately equal: place in center
      newouterposition.legends(i, 2) = newposition.axes(i, 2) + newposition.axes(i, 4)/2 - newouterposition.legends(i, 4)/2;
   end
   
   % finish position (no need for tightinset)
   newposition.legends(i, 1) = newouterposition.legends(i, 1);
   newposition.legends(i, 2) = newouterposition.legends(i, 2);
   newposition.legends(i, 3) = newouterposition.legends(i, 3);
   newposition.legends(i, 4) = newouterposition.legends(i, 4);
end



% *******************************************************************************************************
% apply changes

% axes
setall(handles.axes, 'OuterPosition', newouterposition.axes);
setall(handles.axes, 'Position',      newposition.axes);

% colorbars
setall(handles.colorbars, 'OuterPosition', newouterposition.colorbars);
setall(handles.colorbars, 'Position',      newposition.colorbars);

% axes
setall(handles.legends, 'OuterPosition', newouterposition.legends);
setall(handles.legends, 'Position',      newposition.legends);

% redraw figure (just to make sure)
refresh(fig); drawnow;
pause(min(internalsettings.maxwaittime, internalsettings.waittime));



% *******************************************************************************************************
% save figure

% save
switch lower(filetype)
   case{'png'}
      print(fig, '-dpng', sprintf('-r%i',internalsettings.resolution), [filename,'.png']);
   case{'pdf'}
      print(fig, '-dpdf', sprintf('-r%i',internalsettings.resolution), [filename,'.pdf']);
   case{'jpg'}
      print(fig, sprintf('-djpeg%.0f', internalsettings.jpegqual), sprintf('-r%i',internalsettings.resolution), [filename,'.jpg']); 
   case{'eps'}
      print(fig, '-depsc2', sprintf('-r%i',internalsettings.resolution), [filename,'.eps']);
   case{'epsbw'}
      print(fig, '-deps2', sprintf('-r%i',internalsettings.resolution), [filename,'.eps']);
   case{'+fig'}   
      print(fig, '-dpng', sprintf('-r%i',internalsettings.resolution), [filename,'.png']);
      print(fig, '-depsc2', sprintf('-r%i',internalsettings.resolution), [filename,'.eps']);
      saveas(fig, filename, 'fig');
   otherwise
      print(fig, '-dpng', sprintf('-r%i',internalsettings.resolution), [filename,'.png']);
      print(fig, '-depsc2', sprintf('-r%i',internalsettings.resolution), [filename,'.eps']);
end

% use external conversion tools (produce far better results in some cases)
if internalsettings.conv_png
   system(sprintf(internalsettings.cmd_png2eps, filename, filename)); % convert
   if strcmpi(filetype, 'png') % this means png is just temporary => remove temp file
      [~,~] = system(sprintf(internalsettings.cmd_png2eps_rm, filename));
   end
end

end



% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% get values of FIELD for all object HANDLES (numerical values only; vector length LEN)
function values = getall(handles, field, len)
values = zeros(length(handles), len);
for i = 1: length(handles)
   values(i, :) =  get(handles(i), field);
end
end


% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% set FIELD for all object HANDLES
function setall(handles, field, values)
if ischar(values)
   for i = 1: length(handles)
      set(handles(i), field, values); % string
   end
else
   for i = 1: length(handles)
      set(handles(i), field, values(i,:)); % numerical
   end
end
end


% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% find unique elements of SET with tolerance TOL (quick&dirty, but does not have to be fast)
function set_el = unique_tol(set, tol)
set = set(:);
set_el = [];
for i = 1:length(set)
   % add element to our list
   set_el = [set_el; set(1)]; %#ok<AGROW>
   % find and remove element in set
   set( (set >= set(1)-tol)&(set <= set(1)+tol) ) = [];
   % was this the last element?
   if isempty(set); break; end
end 
end


% *******************************************************************************************************
% *******************************************************************************************************
% *******************************************************************************************************
% make sure settings are backwards compatible (as much as possible)
function internalsettings = backward_compatibility(internalsettings)
% FULLY SUPPORTED
if isfield(internalsettings, 'subplotgap') % gap between two subplots / 2
   warning('savefigure:paramobsolete', 'Field subplotgap is obsolete. Please use axisgap.');
   internalsettings.axisgap(1:end) = internalsettings.subplotgap / 2;
end 
if isfield(internalsettings, 'addspace_x3d')
   warning('savefigure:paramobsolete', 'Field addspace_x3d is obsolete. Please use axisgap3d');
   internalsettings.axisgap3d = internalsettings.axisgap;
   internalsettings.axisgap3d(2) = internalsettings.axisgap3d(2) * internalsettings.addspace_x3d;
end
if isfield(internalsettings, 'colorbarsize')
   warning('savefigure:paramobsolete', 'Field colorbarsize is obsolete. Please use colorbarwidth.');
   internalsettings.colorbarwidth = internalsettings.colorbarsize;
end 
if isfield(internalsettings, 'maxcolbarsize')
   warning('savefigure:paramobsolete', 'Field maxcolbarsize is obsolete. Please use maxcolbarwidth.');
   internalsettings.maxcolbarwidth = internalsettings.maxcolbarsize;
end
% PARTIALLY SUPPORTED
%     can't replicate the "times axes width / axes height" property with the current implementation
if isfield(internalsettings, 'addspace_xlabel')
   warning('savefigure:paramobsolete', 'Field addspace_xlabel is no longer supported; results will differ. Please use axisgap(1).');
   internalsettings.axisgap(1) = internalsettings.axisgap(1) + internalsettings.addspace_xlabel;
end 
if isfield(internalsettings, 'addspace_ylabel')
   warning('savefigure:paramobsolete', 'Field addspace_ylabel is no longer supported; results will differ. Please use axisgap(2).');
   internalsettings.axisgap(2) = internalsettings.axisgap(2) + internalsettings.addspace_ylabel;
end 
if isfield(internalsettings, 'addspace_colorb')
   warning('savefigure:paramobsolete', 'Field addspace_colorb is no longer supported; results will differ. Please use axisgap(3).');
   internalsettings.axisgap(3) = internalsettings.axisgap(3) + internalsettings.addspace_colorb;
end
% NOT SUPPORTED ANY MORE
if isfield(internalsettings, 'smallfig_start')
   warning('savefigure:paramobsolete', 'Field smallfig_start is no longer needed and will be ignored.');
end 
if isfield(internalsettings, 'smallfig_frac')
   warning('savefigure:paramobsolete', 'Field smallfig_frac is no longer needed and will be ignored.');
end 
end



% *******************************************************************************************************
% *******************************************************************************************************
% TEST SCRIPT
% 
% close all
% % XY-PLOTS (legend position approx. subplot position)
% x = linspace(0, 8*pi,1000);
% y = 100 * sin(x);
% %     1 x 1
% figure; plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'NorthEast'); setlabels('TITLE', 'x', 'y');
% savefigure(gcf, 'savefigure_test-xy-1x1', 'png');
% %     2 x 1 
% figure;
% subplot(2,1,1); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(t)'}, 'North'); setlabels('TITLE 1', 'x', 'y');
% subplot(2,1,2); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(t)'}, 'South'); setlabels('TITLE 2', 'x', 'y');
% savefigure(gcf, 'savefigure_test-xy-2x1', 'png');
% %     1 x 2 
% figure;
% subplot(1,2,1); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'West'); setlabels('TITLE 1', 'x', 'y');
% subplot(1,2,2); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'East'); setlabels('TITLE 2', 'x', 'y');
% savefigure(gcf, 'savefigure_test-xy-1x2', 'png');
% %     2 x 2 
% figure;
% subplot(2,2,1); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'NorthWest'); setlabels('TITLE 1', 'x', 'y');
% subplot(2,2,2); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'NorthEast'); setlabels('TITLE 2', 'x', 'y');
% subplot(2,2,3); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'SouthWest'); setlabels('TITLE 3', 'x', 'y');
% subplot(2,2,4); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'SouthEast'); setlabels('TITLE 4', 'x', 'y');
% savefigure(gcf, 'savefigure_test-xy-2x2', 'png');
% %     3 x 2 
% figure;
% subplot(3,2,1); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'NorthWest'); setlabels('TITLE 1', 'x', 'y');
% subplot(3,2,2); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'NorthEast'); setlabels('TITLE 2', 'x', 'y');
% subplot(3,2,3); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'West'); setlabels('TITLE 3', 'x', 'y');
% subplot(3,2,4); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'East'); setlabels('TITLE 4', 'x', 'y');
% subplot(3,2,5); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'SouthWest'); setlabels('TITLE 4', 'x', 'y');
% subplot(3,2,6); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'SouthEast'); setlabels('TITLE 4', 'x', 'y');
% savefigure(gcf, 'savefigure_test-xy-3x2', 'png');
% %     3 x 3 
% figure;
% subplot(3,3,1); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'NorthWest'); setlabels('TITLE 1', 'x', 'y');
% subplot(3,3,2); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'North'); setlabels('TITLE 2', 'x', 'y');
% subplot(3,3,3); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'NorthEast'); setlabels('TITLE 3', 'x', 'y');
% subplot(3,3,4); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'West'); setlabels('TITLE 4', 'x', 'y');
% subplot(3,3,5); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'Best'); setlabels('TITLE 5', 'x', 'y');
% subplot(3,3,6); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'East'); setlabels('TITLE 6', 'x', 'y');
% subplot(3,3,7); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'SouthWest'); setlabels('TITLE 7', 'x', 'y');
% subplot(3,3,8); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'South'); setlabels('TITLE 8', 'x', 'y');
% subplot(3,3,9); plot(x,y); xlim(xyzlimits(x)); ylim(xyzlimits(y));
% setlegend({'sin(x)'}, 'SouthEast'); setlabels('TITLE 9', 'x', 'y');
% savefigure(gcf, 'savefigure_test-xy-3x3', 'png');
% % XYZ-PLOTS
% x = linspace(0, 8*pi,100);
% y = sin(x);
% z = cos(x') * sin(y).^2;
% %     1 x 1
% figure; surf(x,y,z','FaceColor','interp','EdgeColor','none'); xlim(xyzlimits(x)); ylim(xyzlimits(y)); zlim(xyzlimits(z));
% setlegend({'cos(x) * sin(y)^2'}, 'NorthWest'); setcolorbar(); setlabels('TITLE', 'x', 'y', 'z'); view([-10,15]);
% savefigure(gcf, 'savefigure_test-xyz-1x1', 'png');
% %     2 x 2
% figure; 
% subplot(2,2,1); surf(x,y,z','FaceColor','interp','EdgeColor','none'); xlim(xyzlimits(x)); ylim(xyzlimits(y)); zlim(xyzlimits(z));
% setlegend({'cos(x) * sin(y)^2'}, 'NorthWest'); setcolorbar(); setlabels('TITLE 1', 'x', 'y', 'z'); view([-45,30]); 
% subplot(2,2,2); surf(x,y,z','FaceColor','interp','EdgeColor','none'); xlim(xyzlimits(x)); ylim(xyzlimits(y)); zlim(xyzlimits(z));
% setlegend({'cos(x) * sin(y)^2'}, 'NorthEast'); setcolorbar(); setlabels('TITLE 2', 'x', 'y', 'z'); view([-30,30]);
% subplot(2,2,3); surf(x,y,z','FaceColor','interp','EdgeColor','none'); xlim(xyzlimits(x)); ylim(xyzlimits(y)); zlim(xyzlimits(z));
% setlegend({'cos(x) * sin(y)^2'}, 'SouthWest'); setcolorbar(); setlabels('TITLE 3', 'x', 'y', 'z'); view([-15,30]);
% subplot(2,2,4); surf(x,y,z','FaceColor','interp','EdgeColor','none'); xlim(xyzlimits(x)); ylim(xyzlimits(y)); zlim(xyzlimits(z));
% setlegend({'cos(x) * sin(y)^2'}, 'SouthEast'); setcolorbar(); setlabels('TITLE 4', 'x', 'y', 'z'); view([0,30]);
% savefigure(gcf, 'savefigure_test-xyz-2x2', 'png');
% %
% close all