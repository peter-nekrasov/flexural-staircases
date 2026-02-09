zk = 1.2;
d = 1.2;
nu = 0.3;

kappa = 0.3+0.3*1i;


l=2; N = 40; a = 15; M = 1e4;
ns = (0:N).';
sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn = cat(3,sn1,sn2);

ts = linspace(-1,1,501);

% ts = linspace(-10,-1,1000);
% ts = 10.^(ts);

src = []; src.r = [1;0]; src.n = src.r; src.d = [0;1]; src.d2 = [1;0];
src.data = [1;1];
targ = []; targ.r = [cos(ts);sin(ts)]; targ.n = targ.r; targ.d = [-sin(ts);cos(ts)];

ising = 0;
kern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate',kappa,d,sn,[],[],l,ising,nu);
kern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log',nu);

us = kern(src,targ);
u11 = us(1:2:end,1:2:end);
u21 = us(2:2:end,1:2:end);
u12 = us(1:2:end,2:2:end);
u22 = us(2:2:end,2:2:end);

figure(1); t = tiledlayout(2,2);

% nexttile
% plot(log10(ts),(abs(u11))); %,ts,imag(u11));
% 
% nexttile
% plot(log10(ts),(abs(u12))); %,ts,imag(u12));
% 
% nexttile
% plot((ts),(abs(u21))); %,ts,imag(u21));
% 
% nexttile
% plot((ts),(abs(u22))); %,ts,imag(u22));

nexttile
plot(log10(ts),log10(abs(u11))); %,ts,imag(u11));

nexttile
plot(log10(ts),log10(abs(u12))); %,ts,imag(u12));

nexttile
plot(log10(ts),log10(abs(u21))); %,ts,imag(u21));

nexttile
plot(log10(ts),log10(abs(u22))); %,ts,imag(u22));

%%

% kern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth',nu);
% us = kern(src,targ);
% plot(log10(ts),log10(abs(us))); %,ts,imag(u11));
