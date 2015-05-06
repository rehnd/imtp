function [vars, plots] = setup()
%Initialization of variables

% Clean up and add plots to path
clear; close all; clc
addpath('./plots','./prop');

vars.n  = 501;                       % Number of grid points points 
vars.xr = [-5.6, 5.6];               % Range of x values
vars.h  = (vars.xr(2)-vars.xr(1))... % Grid Spacing
            /(vars.n-1);

plots.eigsts = 1;       % Plot initial Eigenstates of H (0 or 1)
plots.showp  = 'on';   % Show plots: 'on' or 'off'
plots.savep  = 1;       % Save plots (0 or 1)
plots.showm  = 'off';   % Show movies: 'on' or 'off'
plots.savem  = 0;       % Save movies (0 or 1)
plots.fn = {'ex_eu_it','ex_eu_rt',...
            'im_eu_it','im_eu_rt',...
            'cr_ni_it','cr_ni_rt',...
            'so_df_it','so_df_rt'};% Cell array of figure names

end