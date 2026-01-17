% chnkr = refine(chnkr,struct('maxchunklen',0.1));
% chnkr = chunkerfunc(@(t) starfish(t,3));
% chnkr = 0.2*chnkr + [d/4;0];
% 
% cgrph = tochunkgraph(chnkr);


zk = 1.2;
nu = 0.3;

nnode = 62;
ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa = ts + amp*1i*sin(ts*d);
xip = 1 + amp*1i*d*cos(ts*d);
ws = ws*xip;

nplot = 240;
nplot = 60;
xx = linspace(-4*d, 4*d,nplot);
yy = linspace(0, 6*d,3*nplot/4) - 1.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
iout = chunkgraphinregion(cgrph,targmod)==1;
src = []; src.r = [[0;-2], [0;2]]; %src.n = [1;0];
    
targout = []; targout.r = targ.r(:,iout);
%%
[sys,sn,l] = clamped_mat(chnkr,zk,nu,kappa,d);
utot = clamped_scatter(src,targout,chnkr,1,0,zk,nu,kappa,d,ws,sys,sn,l);



%%

chnkrs = [];
for i = (-6:6)
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = utot(:,1);

figure(1);clf
h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';
hold on
scatter(src.r(1,1),src.r(2,1),400,'r.')
plot(chnkrs,'k.','markersize',15)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');
% %%
figure(2);clf
us(iout) = utot(:,2);
h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
hold on
scatter(src.r(1,2),src.r(2,2),400,'r.')
plot(chnkrs,'k.','markersize',15)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');



%%
a = 3;
b = 1/a/2;
t0 =-5;
t1 = 5;


f = @(t) flat_interface(t, a, b, t0, t1);
nch = 50;
xrad = -log(eps)/abs(real(zk)) + max([abs(t0),abs(t1)]);
xmin = -xrad;
xmax =  xrad;
cparams = [];
cparams.ta = xmin;
cparams.tb = xmax;
cparams.ifclosed = 0;
chnkr2 = chunkerfuncuni(f, nch, cparams);
% for i = 1:2
% chnkr2 = refine(chnkr2,struct('splitchunks',[nch/2,nch/2+1] )); nch = nch+2;
% chnkr2 = sort(chnkr2);
% end
chnkr2 = [0,-1;1,0]*chnkr2;
%%
% u_tr = clamped_scatter(src,chnkr2,chnkr,1,1,zk,nu,kappa,d,ws,sys,sn,l);
%%
% l=2; N = 40; a = 15; M = 1e4;
% kap = pi/d;
% sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kap,(exp(1i*kap*d)),a,M,l+1);
% 
% ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval_trx',kap,d,sn,s0_l,sn_l,l,1,nu);
% 
% u_tr = ikern(struct('r',[1;0],'n',[1;1]/sqrt(2),'d',[1;-1],'d2',[0.3;0]),chnkr2);
% u = u_tr(4:4:end,1);
% %%
% u = skern_0(src,chnkr2);
% u = u(:,1);


% plot(real(chnkr2.r(2,:)), log10(abs(u)),'.')
% 
% plot(real(chnkr2.r(2,:)), imag(chnkr2.r(2,:)),'.')


%%
src.r = [1;0.];

tic;
u_tr = clamped_scatter(src,chnkr2,chnkr,1,1,zk,nu,kappa,d,ws,sys,sn,l);
toc;
h = 1e-1;
u_tr_p = clamped_scatter(src,chnkr2+[h;0],chnkr,1,1,zk,nu,kappa,d,ws,sys,sn,l);
u_tr_m = clamped_scatter(src,chnkr2+[-h;0],chnkr,1,1,zk,nu,kappa,d,ws,sys,sn,l);

% src.r = [1;0];
% u_tr = direct_rhs(src,chnkr2,zk);
% u_tr_p = direct_rhs(src,chnkr2+[h;0],zk);
% u_tr_m = direct_rhs(src,chnkr2+[-h;0],zk);

%%
a = (u_tr_p -  u_tr_m)/2/h;

[norm(a(1:4:end,:) - u_tr(2:4:end,:)),...
norm(a(2:4:end,:) - u_tr(3:4:end,:)),...
norm(a(3:4:end,:) - u_tr(4:4:end,:))]


i = 2;
plot(real(chnkr2.r(2,:)), log10(abs(a(i:4:end,:) - u_tr(i+1:4:end,:))),'.')

% b = abs(a(i:4:end,:) - u_tr(i+1:4:end,:));
% c = b(abs(chnkr2.r(2,:))< 4).'

b = (u_tr_p +  u_tr_m - 2*u_tr)/h^2;
[norm(b(1:4:end,:) - u_tr(3:4:end,:)),norm(b(2:4:end,:) - u_tr(4:4:end,:))]
% i = 2;
% plot(real(chnkr2.r(2,:)), log10(abs(b(i:4:end,:) - u_tr(i+2:4:end,:))),'.')

% plot(real(chnkr2.r(2,:)), log10(abs(b(i:4:end,:) - a(i+1:4:end,:))),'.')
% figure(3);clf
% plot(real(chnkr2.r(2,:)), real( u_tr(3:4:end,1)),'.')
