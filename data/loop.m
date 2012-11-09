function mpc = loop
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 1;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	2	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	3	1	0	0	0	0	1	1	0	345	1	1.1	0.9;
	4	2	0	0	0	0	1	1	0	345	1	1.1	0.9;
	5	3	0	0	0	0	1	1	0	345	1	1.1	0.9;
	6	1	100	0	0	0	1	1	0	345	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	
mpc.gen = [
	4	0	0	300	-300	1	100	1	100	0	zeros(1,11);
	5	100	0	300	-300	1	100	1	100	0	zeros(1,11);
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	1	4	0	0.10	0	200	250	250	0	0	1	-360	360;
	2	5	0	0.10	0	200	250	250	0	0	1	-360	360;
	3	6	0	0.10    0	200	250	250	0	0	1	-360	360;
	4	5	0	0.10	0	2.5	250	250	0	0	1	-360	360;
	4	6	0	0.10	0	200	250	250	0	0	1	-360	360;
	5	6	0	0.10	0	200	250	250	0	0	1	-360	360;
];


%% generator cost data
%	1	startup	shutdown	n	x1	y1	...	xn	yn
%	2	startup	shutdown	n	c(n-1)	...	c0
mpc.gencost = [
	2	0	0	3	1.5	60	0;
	2	0	0	3	1	40	0;
];
