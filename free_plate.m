zk = 1.2;
d = 1.2;
nu = 0.3;

kappa = pi/d;

nplot = 80;
xx = linspace(-1.5*d, 1.5*d,nplot);
yy = xx;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];


if true
    cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
    nch = 20; A = 1;
    chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
    chnkr = reverse(chnkr);
    wtarg = cos_func(targ.r(1,:),d,A) ;
    iout = targ.r(2,:) > wtarg(2,:);
    src = []; src.r = [0;-2]; src.n = [1;0];
    % src = []; src.r = [0;2]; src.n = [1;0];
else
    chnkr = chunkerfunc(@starfish,struct('eps',1e-10)); chnkr = 0.25*chnkr;
    targmod = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
    iout = ~chunkerinterior(chnkr,targmod);
    src = []; src.r = [0;0]; src.n = [1;0];
end

targout = []; targout.r = targ.r(:,iout);

%%

l=2; N = 40; a = 15; M = 1e4;
ns = (0:N).';
sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn = cat(3,sn1,sn2);

skern = kernel('l','s');
s2trkern = kernel([kernel('l','s');kernel('l','sp')]);

ht = 1.02*d; hb = -1.02*d;
[pxys_l, cs_l] = build_pxys(zk,kappa,d,ht,hb,skern,s2trkern,l,40);

%%


ising = 1;
fkern1 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate',kappa,d,sn,pxys_l,cs_l,l,ising,nu);
double = @(s,t) chnk.lap2dquas.kern(s,t,'d',kappa,d,pxys_l,cs_l,l,ising);
hilbert = @(s,t) chnk.lap2dquas.kern(s,t,'hilb',kappa,d,pxys_l,cs_l,l,ising);

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

start = tic;
sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double, opts);
H = chunkermat(chnkr, hilbert, opts2);     

sysmat = zeros(2*chnkr.npt);
sysmat(1:2:end,1:2:end) = sysmat1(1:4:end,1:2:end) + sysmat1(3:4:end,1:2:end)*H  - 2*((1+nu)/2)^2*D*D;
sysmat(2:2:end,1:2:end) = sysmat1(2:4:end,1:2:end) + sysmat1(4:4:end,1:2:end)*H;
sysmat(1:2:end,2:2:end) = sysmat1(1:4:end,2:2:end) + sysmat1(3:4:end,2:2:end);
sysmat(2:2:end,2:2:end) = sysmat1(2:4:end,2:2:end) + sysmat1(4:4:end,2:2:end);

D = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];  % jump matrix 
D = kron(eye(chnkr.npt), D);

sys =  D + sysmat;
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)


%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa,d,sn,pxys_l,cs_l,l,ising);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa,d,sn,pxys_l,cs_l,l,ising,nu);

rhs = -bskern(src,chnkr);

% Solving linear system
sol = sys\rhs;


ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval',kappa,d,sn,pxys_l,cs_l,l,0,nu);
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

uin = skern(src,targout);
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
quiver(chnkr)
hold on
plot(chnkrs,'k.')
scatter(src.r(1,:),src.r(2,:))
h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';
colorbar
hold off
axis equal



function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end