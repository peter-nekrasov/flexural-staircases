zk = 1.2;
d = 1.2;
nu = 0.3;

kappa = 0.3+0.3*1i;


l=2; N = 40; a = 15; M = 1e4;
ns = (0:N).';
sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn = cat(3,sn1,sn2);

% ts = linspace(-10,-1,1000);
% ts = 10.^(ts);

nplot = 200;
ts = linspace(0,pi/2,nplot);
t0 = ts(nplot/2);

src = []; src.r = [cos(t0);sin(t0)]; src.n = src.r; src.d = [-sin(t0);cos(t0)]; src.d2 = [-cos(t0);-sin(t0)];
src.data = [1;1];
targ = []; targ.r = [cos(ts);sin(ts)]; targ.n = targ.r; targ.d = [-sin(ts);cos(ts)];

ising = 0;
[val, grad, hess, third, fourth,fifth] = chnk.flex2dquas.green(src.r,targ.r,zk,kappa,d,sn,l,0);  
% kern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log',nu);

figure(1); clf;
t = tiledlayout(1,2);

nexttile
plot(ts,real(val)); %,ts,imag(u11));

nexttile
plot(ts,imag(val)); %,ts,imag(u11));

figure(2); clf;
t = tiledlayout(2,2);

nexttile
plot(ts,real(grad(:,:,1))); %,ts,imag(u11));

nexttile
plot(ts,imag(grad(:,:,1))); %,ts,imag(u11));

nexttile
plot(ts,real(grad(:,:,2))); %,ts,imag(u11));

nexttile
plot(ts,imag(grad(:,:,2))); %,ts,imag(u11));

figure(3); clf;
t = tiledlayout(3,2);

nexttile
plot(ts,real(hess(:,:,1))); %,ts,imag(u11));

nexttile
plot(ts,imag(hess(:,:,1))); %,ts,imag(u11));

nexttile
plot(ts,real(hess(:,:,2))); %,ts,imag(u11));

nexttile
plot(ts,imag(hess(:,:,2))); %,ts,imag(u11));

nexttile
plot(ts,real(hess(:,:,3))); %,ts,imag(u11));

nexttile
plot(ts,imag(hess(:,:,3))); %,ts,imag(u11));


% nexttile
% plot(log10(ts),log10(abs(u11))); %,ts,imag(u11));
% 
% nexttile
% plot(log10(ts),log10(abs(u12))); %,ts,imag(u12));
% 
% nexttile
% plot(log10(ts),log10(abs(u21))); %,ts,imag(u21));
% 
% nexttile
% plot(log10(ts),log10(abs(u22))); %,ts,imag(u22));

%%

% kern =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth',nu);
% us = kern(src,targ);
% plot(log10(ts),log10(abs(us))); %,ts,imag(u11));
