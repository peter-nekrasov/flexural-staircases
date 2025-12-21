%%
addpath(genpath('../../flexural-staircases'))
zk = 1;
d = 1.2;
nu = 0.3; 

amp = -0.3;
nleg = 32;
nleg = 64;

xs = cos((2*(1:nleg)-1)/2/nleg*pi);

tmin = zk+0.1; tmax = pi/d;
% tmin = 1.02; tmax = pi/sys.d;
% tmin = 1.001; tmax = 1.01;
tr = (tmax-tmin)*(xs+1)/2+tmin;
kappa = tr;

nkappa = length(kappa);


nplot = 240;
xx = linspace(-6*d, 6*d,nplot);
yy = xx;
yy = linspace(0, 4*d,nplot/2) - 1.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
targmod = targ;
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);

cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
cparams.maxchunklen = 2/zk;cparams.ifclosed = 1;cparams.eps = 1e-6;
nch = 20; A = 1;
chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
chnkr = reverse(chnkr);

wtarg = cos_func(targmod.r(1,:),d,A) ;
iout = targmod.r(2,:) > wtarg(2,:);
targout = []; targout.r = targmod.r(:,iout);
targout_0 = []; targout_0.r = targ.r(:,iout);

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
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

dets = zeros(nkappa,1);

for i = 1:nkappa
    dets(i) = det(2*squeeze(sys(i,:,:)));
end

%%
T = cos((0:(nleg-1)).' .*acos(xs(:).')).';

c_cheb = T\dets;

cs = c_cheb/c_cheb(end);

B = .5*ones(nleg-1,2);
A = spdiags(B,[-1,1],nleg-1,nleg-1);
A(1,2) = 1/sqrt(2);A(2,1) = 1/sqrt(2);
en = zeros(1,nleg-1); en(end)=1;
cs(1) = sqrt(2)*cs(1);

B = A - .5*cs(1:nleg-1)*en;


rts = eig(B);

rts = rts(abs(rts)<1);

rts= rts(abs(imag(rts))<1e-3);

rts = (tmax-tmin)*(rts+1)*.5 + tmin;

figure(5)
plot(kappa,abs(dets))

figure(4);clf
plot(rts,'o')
title('Poles','Interpreter','latex')
set(gca,'fontsize',16)

%%


kappa_rt = real(rts); nkappa = 1;
l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa_rt,(exp(1i*kappa_rt*d)),a,M,l+1);

fkern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate',kappa_rt,d,sn,[],[],l,ising);

curv = signed_curvature(chnkr);
curv = curv(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkern, opts);
sys = reshape(sys,nkappa,2*chnkr.npt,2*chnkr.npt);

sys = sys - reshape(0.5*eye(2*chnkr.npt),1,2*chnkr.npt,2*chnkr.npt);
sys(:,2:2:end,1:2:end) = sys(:,2:2:end,1:2:end) + reshape(curv.*eye(chnkr.npt),1,chnkr.npt,chnkr.npt);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

[u,sig,v] = svd(squeeze(sys));
dens = v(:,end);
sig(end,end)/sig(1,1)

%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa_rt,d,sn,s0_l,sn_l,l,1);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa_rt,d,sn,s0_l,sn_l,l,1,nu);

% Solving linear system
sol = dens;



ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval',kappa,d,sn,[],[],l,0);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');

wts = repmat(chnkr.wts(:).',2,1);
nt = size(targout.r,2);

start1 = tic;
nbatch = ceil(2e5/chnkr.npt);
ntout = size(targout.r,2);

nshiftout = nshift(iout);
for i = 1:ceil(ntout/nbatch)
    iuse = ((i-1)*nbatch+1):min(ntout,i*nbatch);
    targi = []; targi.r = targout.r(:,iuse);
    nti = length(iuse);
    
    gevalmat_0 = chunkerkernevalmat(chnkr,ikern_0,targi,opts);
    gevalmat = ikern(chnkr,targi).* wts(:).';
    
    gevalmat = reshape(gevalmat,nkappa, nti, []);
    gevalmat = gevalmat + reshape(gevalmat_0,1,nti, []);
    gevalmat = exp(1i*kappa(:).*nshiftout(iuse)*d) .* gevalmat;
    % gevalmat = exp(1i*kappa(:).*nshift(iout)*d) .* gevalmat;
    
    gevalmat = reshape(permute(gevalmat, [2,1,3]), nti,[]);
    uscat(iuse,:) = gevalmat*sol;
end
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)



%%
chnkrs = [];
for i = -10:10
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

figure(2);clf
subplot(2,1,1)
us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = uscat;
h = pcolor(X,Y, reshape((imag(us)),size(X))); h.EdgeColor = 'None';
hold on
plot(chnkrs,'k.','markersize',15)
c = colorbar;
c.Label.String = '$\Im v_\xi$';
c.Label.Interpreter = 'latex';
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');

subplot(2,1,2)
h = pcolor(X,Y, reshape((abs(us)),size(X))); h.EdgeColor = 'None';
hold on
plot(chnkrs,'k.','markersize',15)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
xlabel('$x_1$','Interpreter','latex')
ylabel('$x_2$','Interpreter','latex')
c.Label.String = '$|v_\xi|$';
c.Label.Interpreter = 'latex';
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');

exportgraphics(gcf,'clamped_mode.pdf','resolution',200)

function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end