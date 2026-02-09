zk = 1.2;
d = 1.2;
nu = 0.3;

kappa = 0.3+0.3*1i;

nplot = 40;
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

chnkr = makedatarows(chnkr,2);
curv = signed_curvature(chnkr);
kp = arclengthder(chnkr,curv);
kpp = arclengthder(chnkr,kp);

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;

targout = []; targout.r = targ.r(:,iout);

%%

l=2; N = 40; a = 15; M = 1e4;
ns = (0:N).';
sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
sn = cat(3,sn1,sn2);

skern = kernel('l','s');
s2trkern = kernel([kernel('l','s');kernel('l','sp')]);

% ht = 1.02*d; hb = -1.02*d;
% [pxys_l, cs_l] = build_pxys(zk,kappa,d,ht,hb,skern,s2trkern,l,40);

%%
ising = 0;
fkern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate',kappa,d,sn,[],[],l,ising,nu);
fkern_0l =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log',nu);           % build the desired kernel
fkern_0s =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth',nu);           % build the desired kernel

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';

start = tic;
M = chunkermat(chnkr,fkern_0l, opts);
M2 = chunkermat(chnkr,fkern_0s, opts2);
M3 = chunkermat(chnkr,fkern,opts2);
% M3 = 0;

c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + M2 + c0.*curv(:).^2.*eye(chnkr.npt) - 0*eye(chnkr.npt); % extra term shows up for the general problem
M = M + M3 - 0.5*eye(2*chnkr.npt);

sys =  M;
toc(start);
%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa,d,sn,[],[],l,ising);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_bcs',kappa,d,sn,[],[],l,ising,nu);

skern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 's',nu);
bskern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_bcs',nu);

rhs = bskern_0(src,chnkr)+bskern(src,chnkr);

% Solving linear system
sol = sys\rhs;


ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_eval',kappa,d,sn,[],[],l,ising,nu);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval',nu);

start1 = tic;
uscat = chunkerkerneval(chnkr, ikern_0,sol, targout);
toc(start1)
wts = repmat(chnkr.wts(:).',2,1);
uscat = uscat + ikern(chnkr,targout) * (sol .* wts(:));
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

uin = skern_0(src,targout) +skern(src,targout);
utot = uscat(:)-uin(:);



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