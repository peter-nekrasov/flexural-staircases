
ifacc = 1;
d_l = 2;
d_r = 2;


chnkr_l = chunkerfunc(@(t) starfish(t,3,0));
chnkr_l = 0.5*chnkr_l;
% chnkr_l = chnkr_l+ [-d_l/2;0.8];
chnkr_l = merge([chnkr_l+ [-d_l/2;0.8],chnkr_l+ [-d_l/2;-0.8]]);

cgrph_l = tochunkgraph(chnkr_l);


chnkr_r = chunkerfunc(@(t) starfish(t,5,0));
chnkr_r = 0.5*chnkr_r;
chnkr_r = chnkr_r+ [-d_r/2;-0.8];
% chnkr_r = merge([chnkr_r+ [-d_r/2;0.8],chnkr_r+ [-d_r/2;-0.8]]);

cgrph_r = tochunkgraph(chnkr_r);


zk = 1.2;
zk = 7;

nu_l = 0.3;
nu_r = nu_l;

% %%
% d_r = d_l;
% chnkr_r = chnkr_l;
% cgrph_r = cgrph_l;
% nu_r = nu_l;

% d_l = d_r;
% chnkr_l = chnkr_r;
% cgrph_l = cgrph_r;
% nu_l = nu_r;


chnkrs_L = [];
chnkrs_R = [];
for i = (-6:6)
    chnkrs_L = [chnkrs_L, chnkr_l.r(:,:) + [i*d_l;0]];
    chnkrs_R = [chnkrs_R, chnkr_r.r(:,:) + [i*d_r;0]];
end
chnkrs_L = chnkrs_L(:,chnkrs_L(1,:)<0);
chnkrs_R = chnkrs_R(:,chnkrs_R(1,:)>0);


%%
a = 3;
b = 1/a/2;
t0 =-5;
t1 = 5;

t0 = -10;
t1 = 10;

f = @(t) flat_interface(t, a, b, t0, t1);
nch = 50;
% nch = 80;
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
% src = []; src.r = [-1.6;-0.]; 
src = []; src.r = [d_r;-0.]; 
src = []; src.r = [0.5;0.5]; 


nplot = 240;
nplot = 60;nplot = 120;
xx = linspace(-4*d_l, 4*d_l,nplot);
yy = linspace(-3*d_l, 3*d_l,3*nplot/4);
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

% targmod = real([mod(targ.r(1,:)+d_l/2,d_l)-d_l/2;targ.r(2,:)]);
targmod = real([mod(targ.r(1,:),d_l)-d_l;targ.r(2,:)]);
iout = chunkgraphinregion(cgrph_l,targmod)==1;
il = iout & (X(:)<0);
targl = []; targl.r = targ.r(:,il);

% targmod = real([mod(targ.r(1,:)+d_r/2,d_r)-d_r/2;targ.r(2,:)]);
targmod = real([mod(targ.r(1,:),d_r)-d_r;targ.r(2,:)]);
iout = chunkgraphinregion(cgrph_r,targmod)==1;
ir = iout & (X(:)>0);
targr = []; targr.r = targ.r(:,ir);
%%
t1 = tic;
[sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l] = free_mat(chnkr_l,zk,nu_l,kappa_l,d_l);
[sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r] = free_mat(chnkr_r,zk,nu_r,kappa_r,d_r);

%%
[rhsmat_l,layermat_l] = precom_free_layer(chnkr_tr,chnkr_tr,chnkr_l,1,zk,nu_l,kappa_l,d_l,sn_l,l_l,s0_l_l,sn_l_l);
[rhsmat_r,layermat_r] = precom_free_layer(chnkr_tr,chnkr_tr,chnkr_r,1,zk,nu_r,kappa_r,d_r,sn_r,l_r,s0_l_r,sn_l_r);
tpre = toc(t1)




%%
t1 = tic;
data = free_scatter(src,chnkr_tr,chnkr_r,1,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r) ...
    - free_scatter(src,chnkr_tr,chnkr_l,1,1,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l);
% data = free_scatter(src,chnkr_tr,chnkr_r,1,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r);
trhs = toc(t1)

%%

dvals = (-1).^(1:4*chnkr_tr.npt);
Amat = diag(dvals);

apply_sys = @(dens) Amat*dens ...
    + free_layer_fast(chnkr_tr,dens,chnkr_tr,chnkr_r,0,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r,rhsmat_r,layermat_r) ...
    - free_layer_fast(chnkr_tr,dens,chnkr_tr,chnkr_l,0,1,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l,rhsmat_l,layermat_l);

% 

t1 = tic;
dens = gmres(apply_sys,data,[],1e-10,100);
tsolve = toc(t1)

%%
t1 = tic;

ul = free_layer_fast(chnkr_tr,dens,targl,chnkr_l,1,0,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l,rhsmat_l);
ur = free_layer_fast(chnkr_tr,dens,targr,chnkr_r,1,0,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r,rhsmat_r);

tlayer = toc(t1)

t1 = tic;
uinl = free_scatter(src,targl,chnkr_l,1,0,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l);
uinr = free_scatter(src,targr,chnkr_r,1,0,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r);
tin = toc(t1)
%%

uscat = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
uscat(il) = ul;
uscat(ir) = ur;
uin = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
uin(il) = uinl;
uin(ir) = uinr;
%%
us = -uscat+uin;

% us = uscat;
figure(2);clf
h = pcolor(X,Y, reshape((imag(us)),size(X))); h.EdgeColor = 'None';
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


skern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 's');
u0 = skern_0(src,targ);

% us = us-u0(:).';
% figure(3);clf
% h = pcolor(X,Y, reshape((imag(us)),size(X))); h.EdgeColor = 'None';
% hold on
% scatter(src.r(1,1),src.r(2,1),400,'r.')
% plot(chnkrs_L(1,:),chnkrs_L(2,:),'k.','markersize',15)
% plot(chnkrs_R(1,:),chnkrs_R(2,:),'k.','markersize',15)
% c = colorbar;
% hold off
% axis equal
% xlim([min(X(:)),max(X(:))])
% ylim([min(Y(:)),max(Y(:))])
% set(gca,'FontSize',18)
% set(gca,'TickLabelInterpreter','latex');
% set(c,'TickLabelInterpreter','latex');

%%
if ifacc
    src = []; src.r = [-d_r;-0.]; 
    t1 = tic;
    % data = free_scatter(src,chnkr_tr,chnkr_r,1,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r) ...
    %     - free_scatter(src,chnkr_tr,chnkr_l,1,1,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l);
    data = free_scatter(src,chnkr_tr,chnkr_r,1,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r);
    trhs = toc(t1)
    
    %%
    
    dvals = (-1).^(1:4*chnkr_tr.npt);
    Amat = diag(dvals);
    
    apply_sys = @(dens) Amat*dens ...
        + free_layer_fast(chnkr_tr,dens,chnkr_tr,chnkr_r,0,1,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r,rhsmat_r,layermat_r) ...
        - free_layer_fast(chnkr_tr,dens,chnkr_tr,chnkr_l,0,1,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l,rhsmat_l,layermat_l);
    
    % 
    
    t1 = tic;
    dens = gmres(apply_sys,data,[],1e-10,100);
    tsolve = toc(t1)
    
    %%
    t1 = tic;
    
    ul = free_layer_fast(chnkr_tr,dens,targl,chnkr_l,1,0,zk,nu_l,kappa_l,d_l,ws_l,sys_l,sn_l,l_l,H_l,s0_l_l,sn_l_l,rhsmat_l);
    ur = free_layer_fast(chnkr_tr,dens,targr,chnkr_r,1,0,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r,rhsmat_r);
    
    tlayer = toc(t1)
    %%
    t1 = tic;
    uinr = free_scatter(src,targr,chnkr_r,1,0,zk,nu_r,kappa_r,d_r,ws_r,sys_r,sn_r,l_r,H_r,s0_l_r,sn_l_r);
    tin = toc(t1)
    %%
    
    
    % us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
    % us(il) = ul;
    % us(ir) = ur;
    % 
    % figure(2);clf
    % h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
    % hold on
    % scatter(src.r(1,1),src.r(2,1),400,'r.')
    % plot(chnkrs_L(1,:),chnkrs_L(2,:),'k.','markersize',15)
    % plot(chnkrs_R(1,:),chnkrs_R(2,:),'k.','markersize',15)
    % c = colorbar;
    % hold off
    % axis equal
    % xlim([min(X(:)),max(X(:))])
    % ylim([min(Y(:)),max(Y(:))])
    % set(gca,'FontSize',18)
    % set(gca,'TickLabelInterpreter','latex');
    % set(c,'TickLabelInterpreter','latex');
    
    us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
    
    us(il) = 0;
    us(ir) = uinr;
    
    figure(10);clf
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
    
    
    us(il) = ul;
    us(ir) = ur-uinr;
    figure(11);clf
    h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';
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

end