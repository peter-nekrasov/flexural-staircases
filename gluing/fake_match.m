

zk = 1.2;


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
chnkr_tr = chunkerfuncuni(f, nch, cparams);
for i = 1:2
chnkr_tr = refine(chnkr_tr,struct('splitchunks',[nch/2,nch/2+1] )); nch = nch+2;
chnkr_tr = sort(chnkr_tr);
end
chnkr_tr = [0,-1;1,0]*chnkr_tr;



src = []; src.r = [2;0]; 

nplot = 240;
nplot = 60;
xx = linspace(-4, 4,nplot);
yy = linspace(-3, 3,3*nplot/4);
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];
il =(X(:)<0);
targl = []; targl.r = targ.r(:,il);

ir = (X(:)>0);
targr = []; targr.r = targ.r(:,ir);
%%


rhs = direct_rhs(src,chnkr_tr,zk);


%%


dvals = (-1).^(1:4*chnkr_tr.npt);
Amat = diag(dvals);


dens = Amat\rhs;

%%
t1 = tic;

% wts = repmat(chnkr_tr.wts(:).',4,1);

% ul = direct_layer(chnkr_tr,targl,zk)*(dens.*wts(:));
% ur = direct_layer(chnkr_tr,targr,zk)*(dens.*wts(:));
% 
fkern = @(s,t) direct_layer(s,t,zk);
% opts = []; opts.forcesmooth = 'false';
% ul = chunkerkerneval(chnkr_tr,fkern,dens,targl,[]);
% ur = chunkerkerneval(chnkr_tr,fkern,dens,targr);

evalmat = chunkerkernevalmat(chnkr_tr,fkern,targl);
ul = evalmat*dens;
evalmat = chunkerkernevalmat(chnkr_tr,fkern,targr);
ur = evalmat*dens;



tlayer = toc(t1)
%%

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
% us(il) = ul;
us(il) = ul ;
us(ir) = ur;
% us(ir) = 0;


figure(2);clf
% us(iout) = us(:,2);
% h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';
hold on
scatter(src.r(1,1),src.r(2,1),400,'r.')
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
src2 = []; src2.r = src.r;
us(il) = chnk.flex2d.kern(zk,src2,targl,'s');
us(ir) = 0;
figure(4);clf
h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
% h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';
c = colorbar;


us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(il) = ul +chnk.flex2d.kern(zk,src2,targl,'s');

us(ir) = ur;
figure(5);clf
% h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';
c = colorbar;

% figure(3);clf
% 
% plot(real(chnkr_tr.r(2,:)), imag(chnkr_tr.r(2,:)),'.')


