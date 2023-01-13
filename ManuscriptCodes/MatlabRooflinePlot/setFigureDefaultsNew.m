%===============================================================================
%===============================================================================
%
% This function sets the default figure settings for the current session of
% Matlab to some attractive plot settings.
%
%===============================================================================
%
% Input:
%   fig (optional): if provided, it will set the default settings for a
%       particular figure. Otherwise, it will set for all new figures in the
%       current session
% 
% Outputs:
%   figSettings: a struct containing the default figure settings that are
%                set by this function
%         color: 7 x 3 array containing the color order for the default
%                colormap
%        rgbcmy: struct containing some good-looking basic colors 
%                (red, green, blue, cyan, magenta, yellow)
%       xxxList: cell arrays containing a list of markers/lines/dashed lines
%                to cycle through
% Note: a = get(groot, 'default') and b = get(groot, 'factory') print a
% struct of the default or factory settings
% 
%===============================================================================
%===============================================================================

function [figSettings, color, rgbcmy, mList, ...
          lineList, dashList, subplot_pos] = setFigureDefaults(varargin)

p = inputParser;

addParameter(p, 'fig_h', groot);
addParameter(p, 'do_close_all', true);
addParameter(p, 'numcolors', -1);

parse(p, varargin{:});
fig = p.Results.fig_h;
ncolors = p.Results.numcolors;

if p.Results.do_close_all
    % must close all open figures, since settings will only apply to new
    % figures
    close all
end

% if nargin == 0
% %     set defaults to graphics root (i.e., all figures for current session)
%     fig = groot;
% elseif nargin > 1
%     error('SetFigureDefaults may only be called with 1 or 0 arguments')
% end

if ncolors > 0
    color = distinguishable_colors(ncolors);
else
    % get standard line color order
    color = lines(7);
end

% default line width and marker size
lWidth = 3;
mSize = 10;

set(fig, 'defaultFigureUnits', 'normalized');
% vector of where the figure will be positioned by default--would recommend
% changing this based on your screen setup
% [left, bottom, width, height]
% loc = [1921, 1, 200, 100]; % this is for Mike's 2 monitor setup
loc = [1.1, 0, 0.8, 0.8];
% loc = [200, 200, 1500, 900]; % this is for a single monitor
% loc = [50, 50, 1000, 600]; % this is for laptop screen

% this is supposed to keep the figure from being relocated somewhere other
% than 'defaultFigurePosition'
% set(gcf,'Resize','off');

% font size for legend and axes labels
fSize = 34;

% text interpreter to enable latex syntax on labels and axes
textInterp = 'latex';

set(fig, 'defaultLineLineWidth', lWidth);
set(fig, 'defaultLineMarkerSize', mSize);

set(fig, 'defaultScatterLineWidth', lWidth);

set(fig, 'defaultFigurePosition', loc);
set(fig, 'defaultAxesFontSize', fSize);
set(fig, 'defaultAxesTickLabelInterpreter', textInterp);
set(fig, 'defaultTextInterpreter', textInterp);

set(fig, 'defaultLegendInterpreter', textInterp);
set(fig, 'defaultLegendFontSize', fSize);

set(fig, 'defaultTextFontSize', fSize);
set(fig, 'defaultTextFontSize', fSize);

% scatter(nan, nan, '+', 'markerEdgeColor', 'k');
% [~, icons] = legend('a');
%
% % marker size in legend
% set(icons(2), 'defaultGroupMarkerSize', mSize);
% set(icons(2), 'defaultGroupLineWidth', lWidth);

figSettings = get(groot, 'default');

rgbcmy         = struct;
rgbcmy.red     = [254 0 0] / 255;
rgbcmy.green   = [0 179 0] / 255;
rgbcmy.blue    = [0 76 232] / 255;
rgbcmy.cyan    = [0 215 255] / 255;
rgbcmy.magenta = [212 38 255] / 255;
rgbcmy.yellow  = [255 167 0] / 255;

mList = {'+', 'o', '*', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h', '.'};
lineList = {'-+', '-o', '-*', '-x', '-s', '-d', '-^', '-v', '->', '-<', '-p',...
            '-h', '-.'};
dashList = {'--+', '--o', '--*', '--x', '--s', '--d', '--^', '--v', '-->',...
            '--<', '--p', '--h', '--.'};

% this contains good positions for gridded subplots
subplot_pos = struct;

posCell_2x2    = cell(4, 1);
posCell_2x2{1} = [0.03, 0.53, 0.30, 0.42];
posCell_2x2{2} = [0.53, 0.53, 0.30, 0.42];
posCell_2x2{3} = [0.03, 0.03, 0.30, 0.42];
posCell_2x2{4} = [0.53, 0.03, 0.30, 0.42];
%     [left bottom width height]
subplot_pos.dim2x2 = posCell_2x2;

posCell_2x2_square    = cell(4, 1);
posCell_2x2_square{1} = [0.11, 0.56, 0.4, 0.4];
posCell_2x2_square{2} = [0.52, 0.56, 0.4, 0.4];
posCell_2x2_square{3} = [0.11, 0.11, 0.4, 0.4];
posCell_2x2_square{4} = [0.52, 0.11, 0.4, 0.4];
%     [left bottom width height]
subplot_pos.dim2x2_square = posCell_2x2_square;
subplot_pos.dim2x2_outerPosition = [1.1 0.1 0.4 0.71];


posCell_3x1    = cell(3, 1);
posCell_3x1{3} = [0.07, 0.08, 0.88, 0.2];
posCell_3x1{2} = [0.07, 0.40, 0.88, 0.2];
posCell_3x1{1} = [0.07, 0.72, 0.88, 0.2];
%     [left bottom width height]
subplot_pos.dim3x1 = posCell_3x1;

posCell_1x3    = cell(3, 1);
posCell_1x3{1} = [0.03, 0.1, 0.30, 0.88];
posCell_1x3{2} = [0.36, 0.1, 0.30, 0.88];
posCell_1x3{3} = [0.69, 0.1, 0.30, 0.88];
%     [left bottom width height]
subplot_pos.dim1x3 = posCell_1x3;

posCell_3x2    = cell(6, 1);
posCell_3x2{1} = [0.03, 0.53, 0.30, 0.42];
posCell_3x2{2} = [0.34, 0.53, 0.30, 0.42];
posCell_3x2{3} = [0.67, 0.53, 0.30, 0.42];
posCell_3x2{4} = [0.03, 0.05, 0.30, 0.42];
posCell_3x2{5} = [0.34, 0.05, 0.30, 0.42];
posCell_3x2{6} = [0.67, 0.05, 0.30, 0.42];
%     [left bottom width height]
subplot_pos.dim3x2 = posCell_3x2;

end

