%%
% addpath(genpath('../../flexural-staircases'))
% zk = 2;
% zk = 1.2;
zk = 3.6;
zk = 7;
nnode = 61;
ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa = ts + amp*1i*sin(ts*d);
xip = 1 + amp*1i*d*cos(ts*d);
ws = ws*xip;


nkappa = length(kappa);


%%

nplot = 240;
% nplot = 60;
xx = linspace(-4*d, 4*d,nplot);
yy = xx;
yy = linspace(0, 6*d,3*nplot/3) - 1.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
% targmod = targ;
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);
%%

src = []; src.r = [[0;-2],[d/2;1.5]];
% src = []; src.r = [[0;-2],[d/2;1]];
% src.r = [0;-2];

chnkrs = [];
for i = (-6:6)
    chnkrs = [chnkrs, chnkrplot + [i*d;0]];
end
chnkrs = merge(chnkrs);
% iout = ~chunkerinterior(chnkrs,targ);
iout = chunkgraphinregion(cgrph,targmod)==1;

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
toc(start)
%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa,d,sn,[],[],l,ising);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'clamped_plate_bcs',kappa,d,sn,[],[],l,ising);

skern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 's');

rhs = -bskern(src,chnkr);

% Solving linear system
sol = 0*rhs;
for i = 1:nkappa
sol(i:nkappa:end,:) = (squeeze(sys(i,:,:))\rhs(i:nkappa:end,:))*ws(i);
cond(squeeze(sys(i,:,:)))
end

%%

uin = skern_0(src,targout_0);
if nkappa == 1
    uin = exp(1i*kappa(:).*nshift(iout).'*d).*skern(src,targout_0)*ws(1);
end

uscat = 0*uin;

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


utot = uscat+uin;



%%
% chnkrs = [];
% for i = -4:4
%     chnkrs = [chnkrs, chnkr + [i*d;0]];
% end
% chnkrs = merge(chnkrs);

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = utot(:,1);

f1=figure(1);clf
f1.Position = [1 1 643 441];
C = reshape(log10(abs(us)/norm(uin(:,1),inf)), size(X)); 
h = pcolor(X,Y,C);
set(gca,'Color','w')
h.EdgeColor = 'None'; 
h.FaceColor = 'texturemap'; 
h.AlphaData = ~isnan(C);
h.FaceAlpha = 'texturemap';
hold on
scatter(src.r(1,1),src.r(2,1),400,'r.')
plot(chnkrs,'k-','LineWidth',2.5)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
vv = sort(abs(C(:)));
clim([min(C(:)),-vv(3)])
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');
% exportgraphics(f1,'clamped_acc.pdf','ContentType','vector','Resolution',300);
% %%
f2=figure(2);clf
f2.Position = [1 1 643 441];
us(iout) = utot(:,2);
C = reshape((imag(us)),size(X));
h = pcolor(X,Y,C);
h.EdgeColor = 'None'; 
h.FaceColor = 'texturemap'; 
h.AlphaData = ~isnan(C);
h.FaceAlpha = 'texturemap';
hold on
scatter(src.r(1,2),src.r(2,2),300,'r.')
plot(chnkrs,'k-','LineWidth',2.5)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');
% exportgraphics(f2,'clamped_sol.pdf','ContentType','vector','Resolution',300);


function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end