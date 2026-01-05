zks = [0.533608321252813,  0.913684049343464, 1.38315377048501];

nplot = 200;
% nplot = 60;
xx = linspace(-4*d, 4*d,nplot);
yy = xx;
yy = linspace(0, 4*d,nplot/2) - 1.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
% targmod = targ;
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);

% wtarg = cos_func(targmod.r(1,:),d,A) ;
% iout = targmod.r(2,:) > wtarg(2,:);
iout = chunkgraphinregion(cgrph,targmod)==1;
targout = []; targout.r = targmod.r(:,iout);
targout_0 = []; targout_0.r = targ.r(:,iout);


%%
figure(10);clf; t = tiledlayout('flow'); t.TileSpacing = 'tight';
figure(11);clf; t1 = tiledlayout('flow'); t1.TileSpacing = 'tight';


for j = 1:length(zks)

zk = real(zks); kappa_rt = kappas(end); nkappa = 1;
if nkappa == 0, return, end
start = tic;
% fkern1 =  @(s,t) chnk.helm2dquas.kern(zk, s, t, 'sprime',quas_param);
fkern1 = kernel('hq','sp',zk,kappa_rt,d);
sysmat1 = chunkermat(chnkr,fkern1);

sys = -0.5*eye(size(sysmat1)) + sysmat1;
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

[u,sig,v] = svd(squeeze(sys));
dens = v(:,end);
sig(end,end)/sig(1,1)
%%

ikern =  kernel('hq','s',zk,kappa_rt,d,1);
% Solving linear system
sol = dens;
ikern_0 = 0*kernel('h','s',zk);

wts = chnkr.wts(:).';

start1 = tic;
uscat = chunkerkerneval(chnkr, ikern_0,dens,targout);
uscat = uscat + ikern.eval(chnkr,targout) * (dens .* wts(:));
uscat = uscat.*exp(1i*kappa_rt*d*nshift(iout).');
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


%%
chnkrs = [];
for i = -10:10
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

figure(10);nexttile()
us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = uscat;
h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
hold on
plot(chnkrs,'k.','markersize',10)
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

figure(11);nexttile()
h = pcolor(X,Y, reshape((abs(us)),size(X))); h.EdgeColor = 'None';
hold on
plot(chnkrs,'k.','markersize',10)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
c.Label.String = '$|v_\xi|$';
c.Label.Interpreter = 'latex';
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');

end
