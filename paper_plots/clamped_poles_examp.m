%%
addpath(genpath('../../flexural-staircases'))
load('clamped_fake_poles1.mat')
% d = 1.2;
nu = 0.3; 


nplot = 240;
nplot = 120;
xx = linspace(-4*d, 4*d,nplot);
yy = xx;
yy = linspace(0, 4*d,nplot/2) - 1.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
% targmod = targ;
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);

% cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
% cparams.maxchunklen = 2/zk;cparams.ifclosed = 1;cparams.eps = 1e-6;
% nch = 20; A = 1;
% chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
% chnkr = reverse(chnkr);
iout = chunkgraphinregion(cgrph,targmod)==1;
targout = []; targout.r = targmod.r(:,iout);
targout_0 = []; targout_0.r = targ.r(:,iout);

%%
figure(2);clf;
t = tiledlayout('flow'); t.Padding = 'loose';
for j = 1:length(zks)
    % for j = 7
    zk = zks(j);
if isempty(poles{j}), continue, end
kappa_rt = real(poles{j}(1)); nkappa = 1;
l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa_rt,(exp(1i*kappa_rt*d)),a,M,l+1);

fkern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate',kappa_rt,d,sn,[],[],l,ising);

curv = signed_curvature(chnkr);
curv = curv(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkfrn, opts);
sys = reshape(sys,nkappa,2*chnkr.npt,2*chnkr.npt);

sys = sys - reshape(0.5*eye(2*chnkr.npt),1,2*chnkr.npt,2*chnkr.npt);
sys(:,2:2:end,1:2:end) = sys(:,2:2:end,1:2:end) + reshape(curv.*eye(chnkr.npt),1,chnkr.npt,chnkr.npt);
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

[u,sig,v] = svd(squeeze(sys));
dens = v(:,end);
sig(end,end)/sig(1,1)

%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa_rt,d,sn,[],[],l,1);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa_rt,d,sn,[],[],l,1,nu);

% Solving linear system
sol = dens;



ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_eval',kappa_rt,d,sn,[],[],l,0);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'clamped_plate_eval');

wts = repmat(chnkr.wts(:).',2,1);
nt = size(targout.r,2);

start1 = tic;
nbatch = ceil(2e5/chnkr.npt);
ntout = size(targout.r,2);
uscat = zeros(ntout,1);
nshiftout = nshift(iout);
for i = 1:ceil(ntout/nbatch)
    iuse = ((i-1)*nbatch+1):min(ntout,i*nbatch);
    targi = []; targi.r = targout.r(:,iuse);
    nti = length(iuse);
    
    gevalmat_0 = chunkerkernevalmat(chnkr,ikern_0,targi,opts);
    gevalmat = ikern(chnkr,targi).* wts(:).';
    
    gevalmat = reshape(gevalmat,nkappa, nti, []);
    gevalmat = gevalmat + reshape(gevalmat_0,1,nti, []);
    gevalmat = exp(1i*kappa_rt(:).*nshiftout(iuse)*d) .* gevalmat;
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

nexttile()
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

end
function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end