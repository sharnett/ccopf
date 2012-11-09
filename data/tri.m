 function mpc = case9
%CASE9    Power flow data for triangle
%   Please see CASEFORMAT for details on the case file format.
%
%   Based on data from Hiskens paper

%   MATPOWER
%   $Id: triangle.m,v 1.11 2012/04/06 17:46:14 by Dano $

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	2	0	0	0	0	1	1.05	64	1	1	1.1	0.9;
	2	2	0	0	0	0	1	1.05	70	1	1	1.1	0.9;
	3	1	100	0	0	0	1	1.05	30	1	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	1	0	300	-300	1	1	1	250	10	0	0	0	0	0	0	0	0	0	0	0;
	2	0.5	1.8	300	-300	1	1	1	300	10	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	2	0	1.0	0	6	250	250	0	0	1	-360	360;
	2	3	0	1.0	0	100	250	250	0	0	1	-360	360;
	1	3	0	1.0	0	100	250	250	0	0	1	-360	360;
];

%%-----  OPF Data  -----%%
%% area data
%	area	refbus
mpc.areas = [
	1	3;
];

%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	1500	0	3	0.11	5	150;
	2	2000	0	3	0.085	1.2	600;
];
