
% This file is taken from the following library: https://github.com/JackDelcaro/MATLAB-graphic-tools
% Check the library documentation for more information

%% DEFAULT GRAPHIC SETTINGS

set(groot,'defaultAxesFontSize',21);                    % Font size for axes
set(groot,'defaultTextFontSize',21);                    % Font size for text
set(groot,'defaultAxesFontName','latex');               % Font name for axes
set(groot,'defaultTextFontName','latex');               % Font name for text
set(groot,'DefaultAxesBox','on');                       % Enable box in graphics
set(groot,'DefaultAxesXGrid','on');                     % Enable grid in graphics
set(groot,'DefaultAxesYGrid','on');                     % Enable grid in graphics
set(groot,'DefaultLineLinewidth',2);                    % Line width for plots
set(groot, 'DefaultStairLineWidth', 2);                 % Line width for stairs

set(0,'DefaultFigureWindowStyle','docked');             % Set figures to docked

set(0, 'defaultAxesTickLabelInterpreter','latex');      % Axes tick label
set(0, 'defaultLegendInterpreter','latex');             % Legend
set(0, 'defaultTextInterpreter','latex');               % Miscellaneous strings
set(0, 'defaultColorBarTickLabelInterpreter', 'latex'); % Color bar ticks
