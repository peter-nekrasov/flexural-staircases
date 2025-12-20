%%
zk = 2;
d = 1.2;

nnode = 61;
ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa = ts + amp*1i*sin(ts*d);
xip = 1 + amp*1i*d*cos(ts*d);
ws = ws*xip;

% kappa = kappa(1);
% ws = ws(1);

nkappa = length(kappa);


%%

nplot = 80;
xx = linspace(-3*d, 3*d,nplot);
yy = xx;
yy = linspace(0, 6*d,nplot) - 1.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];


%%

if true
    cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
    cparams.maxchunklen = 2/zk;cparams.ifclosed = 0;
    nch = 20; A = 1;
    % chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
    chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
    chnkr = reverse(chnkr);
    wtarg = cos_func(targ.r(1,:),d,A) ;
    iout = targ.r(2,:) > wtarg(2,:);
    % src = []; src.r = [0;-3]; src.n = [1;0];
    src = []; src.r = [0;2]; src.n = [1;0];
else
    chnkr = chunkerfunc(@starfish,struct('eps',1e-10)); chnkr = 0.25*chnkr;
    targmod = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
    iout = ~chunkerinterior(chnkr,targmod);
    src = []; src.r = [0;0]; src.n = [1;0];
end

targout = []; targout.r = targ.r(:,iout);

%%

l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);

%%
ising = 1;
fkern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate',kappa,d,sn,[],[],l,ising);

curv = signed_curvature(chnkr);
curv = curv(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkern, opts);
sys = reshape(sys,nkappa,2*chnkr.npt,2*chnkr.npt);

sys = sys - reshape(0.5*eye(2*chnkr.npt),1,2*chnkr.npt,2*chnkr.npt);
sys(:,2:2:end,1:2:end) = sys(:,2:2:end,1:2:end) + reshape(curv.*eye(chnkr.npt),1,chnkr.npt,chnkr.npt);
toc(start)
%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa,d,sn,[],[],l,ising);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_bcs',kappa,d,sn,[],[],l,ising);

skern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 's');

rhs = -bskern(src,chnkr);

% Solving linear system
sol = 0*rhs;
for i = 1:nkappa
sol(i:nkappa:end) = (squeeze(sys(i,:,:))\rhs(i:nkappa:end))*ws(i);
end

%%
ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval',kappa,d,sn,[],[],l,0);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');

wts = repmat(chnkr.wts(:).',2,1);
nt = size(targout.r,2);

start1 = tic;
gevalmat_0 = chunkerkernevalmat(chnkr,ikern_0,targout,opts);
gevalmat = ikern(chnkr,targout).* wts(:).';

gevalmat = reshape(gevalmat,nkappa, nt, []);
gevalmat = gevalmat + reshape(gevalmat_0,1,nt, []);
gevalmat = reshape(permute(gevalmat, [2,1,3]), nt,[]);
uscat = gevalmat*sol;
t2 = toc(start1);

% start1 = tic;
% eval_sing = chunkerkernevalmat(chnkr, ikern_0, targout);
% sol0 = ws*reshape(sol,nkappa,[]);
% uscat = eval_sing * sol0(:);
% toc(start1)
% 
% wts = repmat(chnkr.wts(:).',2,1);
% wts = ws.'.*reshape(wts,1,[]);
% eval_smooth = ikern(chnkr,targout);
% 
% uscat = uscat +  * (sol .* wts(:));
% t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

uin = skern_0(src,targout);
if nkappa == 1
    uin = skern(src,targout);
end
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
% h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';
h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
clim([-0.1,0.1]*0.5)
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