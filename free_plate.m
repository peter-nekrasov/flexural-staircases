zk = 2;
% zk = 0.4;
d = 1.2;
% d = 2*pi;
nu = 0.3;

% kappa = pi/d;
% kappa = 0.5;

nplot = 80;
xx = linspace(-1.5*d, 1.5*d,nplot);
yy = xx;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];


if true
    cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
    cparams = []; cparams.ta = 0; cparams.tb = d;
    nch = 20/2; A = 0.2;
    chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
    chnkr = reverse(chnkr);
    wtarg = cos_func(targ.r(1,:),d,A) ;
    iout = targ.r(2,:) > wtarg(2,:);
    src = []; src.r = [0;-2]; src.n = [1;0];
    coefs = 1;
    src = []; src.r = [[-0.2;-2],[0.2;-2]]; src.n = [[1;0],[1;0]];
    coefs = [1;-1];
    % src = []; src.r = [0;2]; src.n = [1;0];

else
    chnkr = chunkerfunc(@(t) starfish(t,3,0.1),struct('eps',1e-10,'maxchunklen',0.4)); chnkr = 0.25*chnkr;
    chnkr = chnkr*1.3;
    targmod = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
    iout = ~chunkerinterior(chnkr,targmod);
    src = []; src.r = [[-0.01;-0.01], [0.05;0.02]]; %src.n = [1;0];
    coefs = [1;0];
end

targout = []; targout.r = targ.r(:,iout);

%%

l=2; N = 40; a = 15; M = 1e4;
ns = (0:N).';
sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn = cat(3,sn1,sn2);

% skern = kernel('l','s');
% s2trkern = kernel([kernel('l','s');kernel('l','sp')]);
% 
% ht = 1.02*d; hb = -1.02*d;
% % ht = 1.1*d; hb = -1.1*d;
% [pxys_l, cs_l] = build_pxys(0,kappa,d,ht,hb,skern,s2trkern,l,60);

ns = 1:N;
sn_l = 0;
s0_l = 0;
for j = [-l:-1, 1:l]
    sn_l = sn_l + j.^-ns.*exp(1i*kappa*d * j);
    s0_l = s0_l + log(j*d).*exp(1i*kappa*d * j);
end

sn_l = -sn_l + polylog(ns,exp(1i*kappa*d)) + (-1).^-ns.*polylog(ns,exp(-1i*kappa*d));
sn_l = sn_l./ns./d.^ns;
sn_l = sn_l / 2 / pi ;
s0_l = quasi_dual_sum(0,d/2,0,kappa,d).'- chnk.lap2dquas.green([0;0],[0;d/2],kappa,d,0,sn_l,l,1);

%%


ising = 0;
fkern1 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate',kappa,d,sn,s0_l,sn_l,l,ising,nu);
double = @(s,t) chnk.lap2dquas.kern(s,t,'d',kappa,d,s0_l,sn_l,l,ising);
hilbert = @(s,t) chnk.lap2dquas.kern(s,t,'hilb',kappa,d,s0_l,sn_l,l,ising);
opts = [];
opts.sing = 'smooth';

opts2 = [];
opts2.sing = 'smooth';

% building system matrix

start = tic;
sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double, opts);
H = chunkermat(chnkr, hilbert, opts2);     


fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate',nu);
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

sysmat1_0 = chunkermat(chnkr,fkern1, opts);
D_0 = chunkermat(chnkr, double, opts);
H_0 = chunkermat(chnkr, hilbert, opts2); 

% sysmat1 = sysmat1 + sysmat1_0; D = D_0; H = H_0;

sysmat1 = sysmat1 + sysmat1_0; D = D + D_0; H = H + H_0;

% 
% ising = 1;
% fkern1 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate',kappa,d,sn,pxys_l,cs_l,l,ising,nu);
% double = @(s,t) chnk.lap2dquas.kern(s,t,'d',kappa,d,pxys_l,cs_l,l,ising);
% hilbert = @(s,t) chnk.lap2dquas.kern(s,t,'hilb',kappa,d,pxys_l,cs_l,l,ising);
% 
% opts = [];
% opts.sing = 'log';
% 
% opts2 = [];
% opts2.sing = 'pv';
% 
% % building system matrix
% 
% sysmat1 = chunkermat(chnkr,fkern1, opts);
% D = chunkermat(chnkr, double, opts);
% H = chunkermat(chnkr, hilbert, opts2); 

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) + sysmat1(3:4:end,1:2:end)*H  - 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) + sysmat1(4:4:end,1:2:end)*H;
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

% for j = 1:2
%     for k = 1:2
% sysmat(j:2:end,k:2:end) = sysmat(j:2:end,k:2:end) + chnkr.wts(:).';
%     end
% end
D = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];  % jump matrix 
D = kron(eye(chnkr.npt), D);

sys =  D + sysmat;
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)


%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa,d,sn,s0_l,sn_l,l,1);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa,d,sn,s0_l,sn_l,l,1,nu);

rhs = -bskern(src,chnkr)*coefs;

% Solving linear system
sol = sys\rhs;


ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval',kappa,d,sn,s0_l,sn_l,l,0,nu);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval',nu);


dens_comb = zeros(3*chnkr.npt,1);
dens_comb(1:3:end) = sol(1:2:end);
dens_comb(2:3:end) = H*sol(1:2:end);
dens_comb(3:3:end) = sol(2:2:end);

wts = repmat(chnkr.wts(:).',3,1);

start1 = tic;
uscat = chunkerkerneval(chnkr, ikern_0,dens_comb,targout);
uscat = uscat + ikern(chnkr,targout) * (dens_comb .* wts(:));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

uin = skern(src,targout)*coefs;
utot = uscat(:)+uin(:);


%%
chnkrs = [];
for i = -1:1
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = utot;

figure(2);clf
h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';

hold on
plot(chnkrs,'k.')
scatter(src.r(1,:),src.r(2,:))
quiver(chnkr)
% h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
% clim(2e-10*[-1,1])
colorbar
hold off
axis equal

% %%
% figure(3)
% 
% plot(chnkr.r(1,:), sol(1:2:end),'.')
% hold on
% plot(chnkr.r(1,:)+d, sol(1:2:end),'.')
% % plot(chnkr.r(1,:), sol(2:2:end),'.')
% % plot(chnkr.r(1,:)+d, sol(2:2:end),'.')
% 
% hold off

%%
% 
% 
% [xs,ws] = lege.exps(60); 
% xs = 0.5*(xs(:).'+1); ws = 0.5*ws(:).';
% 
% ts = pi/d*[-xs, xs];
% 
% ws = pi/d*[ws,ws]/2/pi;
% 
% xi = ts - 0.3i*sin(ts*d);
% xip = ws.*(1 - 0.3i*d*cos(ts*d));
% 
% skern_0 = kernel('l','s');
% s2trkern = kernel([kernel('l','s');kernel('l','sp')]);
% npxy = 40;
% [pxys, cs] = build_pxys(0,xi,d,ht,hb,skern_0,s2trkern,l,npxy);
% 
% 
% src = []; src.r = [0;0]; src.n = chnkr.n(:,1); src.d = chnkr.d(:,1); src.d2 = chnkr.d2(:,1);
% targ = []; targ.r = [0.4;0.]; targ.n = chnkr.n(:,100); targ.d = chnkr.d(:,100); targ.d2 = chnkr.d2(:,100);
% 
% free_kap = 0;
% 
% 
% vals = zeros(4,2,length(xi));
% for i = 1:length(xi)
% 
%     l=2; N = 40; a = 15; M = 1e4;
% ns = (0:N).';
% sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,xi(i),(exp(1i*xi(i)*d)),a,M,l+1);
% sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,xi(i),(exp(1i*xi(i)*d)),a,M,l+1);
% sn = cat(3,sn1,sn2);
% 
% free_i = chnk.flex2dquas.kern(zk, src, targ, 'free_plate',xi(i),d,sn,pxys_l,cs(:,i),l,1,nu);
% vals(:,:,i) = free_i;
% free_kap = free_kap + free_i * xip(i);
% end
% 
% 
% 
% free_0 = chnk.flex2d.kern(zk, src, targ, 'free_plate',nu);
% 
% err = abs(free_kap - free_0)
% free_kap
% free_0



% figure(1);clf
% i = 1; j = 1;
% plot(ts, [squeeze(real(vals(i,j,:))), squeeze(imag(vals(i,j,:)))],'.')
% % 

% %%
% 
% [U,s,V] = svd(sys);
% figure(1);clf
% plot(log10(diag(s)),'.')
% %%
% v = V(:,end);
% plot(chnkr.r(1,:),imag(v(1:2:end)),'.')
% hold on
% plot(chnkr.r(1,:),imag(v(2:2:end)),'.')
% hold off
% 
% 
% 
% % %%
% a = exp(1i*kappa*chnkr.r(1,:).');
% % a = 1+0.*chnkr.r(1,:).';
% b = H*a;
% 
% c = mean(b./a);
% 
% % norm(b - c*a)
% 
% % plot(chnkr.r(1,:),real(a),'.')
% plot(chnkr.r(1,:),real(b./a),'.')
% hold on
% % plot(chnkr.r(1,:),real(b),'.')
% hold off

function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end