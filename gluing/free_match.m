
d_l = 2;
d_r = 2.5;


chnkr_l = chunkerfunc(@(t) starfish(t,3));
chnkr_l = 0.2*chnkr_l + [-d_l/4;0];

cgrph_l = tochunkgraph(chnkr_l);


chnkr_r = chunkerfunc(@(t) starfish(t,3));
chnkr_r = 0.2*chnkr_r + [d_r/4;-1];

cgrph_r = tochunkgraph(chnkr_r);


zk = 1.2;

nu_l = 0.3;
nu_r = nu_l;

%%
d_r = d_l;
chnkr_r = chnkr_l;
cgrph_r = cgrph_l;
nu_r = nu_l;


%%
a = 3;
b = 1/a/2;
t0 =-5;
t1 = 5;

t0 = -10;
t1 = 10;

f = @(t) flat_interface(t, a, b, t0, t1);
nch = 50;
xrad = -log(eps)/abs(real(zk)) + max([abs(t0),abs(t1)]);
xmin = -xrad;
xmax =  xrad;
cparams = [];
cparams.ta = xmin;
cparams.tb = xmax;
cparams.ifclosed = 0;
chnkr_tr = chunkerfuncuni(f, nch, cparams);
for i = 1:2
    isplit = find(abs(min(chnkr_tr.r(1,:,:))) < 2);
chnkr_tr = refine(chnkr_tr,struct('splitchunks',isplit(:).')); nch = nch+length(isplit);
chnkr_tr = sort(chnkr_tr);
end

% for i = 1:2
%     isplit = find(abs(min(chnkr_tr.r(1,:,:))) < 1);
% chnkr_tr = refine(chnkr_tr,struct('splitchunks',isplit(:).')); nch = nch+length(isplit);
% chnkr_tr = sort(chnkr_tr);
% end

chnkr_tr = [0,-1;1,0]*chnkr_tr;


%%
nnode = 62;
ts = linspace(-pi/d_l,pi/d_l,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa_l = ts + amp*1i*sin(ts*d_l);
xip = 1 + amp*1i*d_l*cos(ts*d_l);
ws_l = ws*xip;

ts = linspace(-pi/d_r,pi/d_r,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa_r = ts + amp*1i*sin(ts*d_r);
xip = 1 + amp*1i*d_r*cos(ts*d_r);
ws_r = ws*xip;

%%
src = []; src.r = [-2;-1]; 

nplot = 240;
nplot = 60;
xx = linspace(-4*d_l, 4*d_l,nplot);
yy = linspace(-3*d_l, 3*d_l,3*nplot/4);
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = real([mod(targ.r(1,:)+d_l/2,d_l)-d_l/2;targ.r(2,:)]);
iout = chunkgraphinregion(cgrph_l,targmod)==1;
il = iout & (X(:)<0);
targl = []; targl.r = targ.r(:,il);

targmod = real([mod(targ.r(1,:)+d_r/2,d_r)-d_r/2;targ.r(2,:)]);
iout = chunkgraphinregion(cgrph_r,targmod)==1;
ir = iout & (X(:)>0);
targr = []; targr.r = targ.r(:,ir);
%%
[sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l] = free_mat(chnkr_l,zk,nu_l,kappa_l,d_l);
[sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r] = free_mat(chnkr_r,zk,nu_r,kappa_r,d_r);


%%
t1 = tic;
% data = free_scatter(src,chnkr_tr,chnkr_r,1,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r) ...
%     - free_scatter(src,chnkr_tr,chnkr_l,1,1,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l);
data = free_scatter(src,chnkr_tr,chnkr_r,1,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r);
trhs = toc(t1)

%%

dvals = (-1).^(1:4*chnkr_tr.npt);
Amat = diag(dvals);


apply_sys = @(dens) Amat*dens ...
    + free_layer(chnkr_tr,dens,chnkr_tr,chnkr_r,0,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r) ...
    - free_layer(chnkr_tr,dens,chnkr_tr,chnkr_l,0,1,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l);

% 
% 
% t1 = tic;
% dens = gmres(apply_sys,data,[],1e-6,100);
% tsolve = toc(t1);

dens = Amat\data;
% dens = data;

%%
t1 = tic;

ul = free_layer(chnkr_tr,dens,targl,chnkr_l,1,0,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l);
ur = free_layer(chnkr_tr,dens,targr,chnkr_r,1,0,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r);

tlayer = toc(t1)
%%
t1 = tic;
uinr = free_scatter(src,targr,chnkr_r,1,0,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r);
tin = toc(t1)
%%

chnkrs_L = [];
chnkrs_R = [];
for i = (-6:6)
    chnkrs_L = [chnkrs_L, chnkr_l.r(:,:) + [i*d_l;0]];
    chnkrs_R = [chnkrs_R, chnkr_r.r(:,:) + [i*d_r;0]];
end
chnkrs_L = chnkrs_L(:,chnkrs_L(1,:)<0);
chnkrs_R = chnkrs_R(:,chnkrs_R(1,:)>0);


us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(il) = ul;
us(ir) = ur;

figure(2);clf
h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
hold on
scatter(src.r(1,1),src.r(2,1),400,'r.')
plot(chnkrs_L(1,:),chnkrs_L(2,:),'k.','markersize',15)
plot(chnkrs_R(1,:),chnkrs_R(2,:),'k.','markersize',15)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));

us(il) = 0;
us(ir) = uinr;

figure(3);clf
h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
hold on
scatter(src.r(1,1),src.r(2,1),400,'r.')
plot(chnkrs_L(1,:),chnkrs_L(2,:),'k.','markersize',15)
plot(chnkrs_R(1,:),chnkrs_R(2,:),'k.','markersize',15)

c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');



% 
% %%
% 
% %%
% u_tr = free_scatter(src,chnkr_tr,chnkr,1,zk,nu,kappa,d,ws,sys,sn,l,H,s0_l,sn_l);
% %%
% % l=2; N = 40; a = 15; M = 1e4;
% % kap = pi/d;
% % sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kap,(exp(1i*kap*d)),a,M,l+1);
% % 
% % ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval_trx',kap,d,sn,s0_l,sn_l,l,1,nu);
% % 
% % u_tr = ikern(struct('r',[1;0],'n',[1;1]/sqrt(2),'d',[1;-1],'d2',[0.3;0]),chnkr2);
% u = u_tr(4:4:end,1);
% % %%
% % u = skern_0(src,chnkr2);
% % u = u(:,1);
% 
% 
%%
figure(5);clf
% plot(real(chnkr_tr.r(2,:)), real(dens(4:4:end)),'.')
plot(real(chnkr_tr.r(2,:)), log10(abs(dens(4:4:end))),'.')
% plot(real(chnkr2.r(2,:)), log10(abs(u)),'.')
% 
% plot(real(chnkr_tr.r(2,:)), imag(chnkr_tr.r(2,:)),'.')


a = reshape(dens,4,16,[]);

[~,~,val2coef] = lege.exps(16);


b = val2coef*squeeze(a(4,:,:));

figure(6);
plot(1:16,log10(b),'-')
    % plot(1:16,log10(abs(b(:,22:22))),'-')