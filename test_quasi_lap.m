zk = 0;
kappa = pi/d-0.1+0.4i;

d = 0.7;

n = 100;
npxy = 40;
thetas = 2*pi*(1:npxy)/npxy;
pxy = 5*d/2 * [cos(thetas);sin(thetas)];
pxys = []; pxys.r = pxy; pxys.n = [cos(thetas);sin(thetas)];

xs = linspace(-d/2,d/2,n);
ys = linspace(-3*d/2,3*d/2,2*n);


rup = [xs(:).'; 0*xs(:).' + 3*d/2];
rdw = [xs(:).'; 0*xs(:).' - 3*d/2];

rlf = [0*ys(:).' - d/2; ys(:).'];
rrt = [0*ys(:).' + d/2; ys(:).'];

rwall = [];
rwall.r = [[0*ys(:).' - d/2; ys(:).'],[0*ys(:).' + d/2; ys(:).']];
rwall.n = [1;0] + 0*rwall.r;


% figure(1); clf
% scatter(rup(1,:), rup(2,:))
% hold on
% scatter(rdw(1,:), rdw(2,:))
% scatter(rlf(1,:), rlf(2,:))
% scatter(rrt(1,:), rrt(2,:))
% scatter(pxy(1,:), pxy(2,:))
% hold off
% axis equal


skern = kernel('l','s');
spkern = kernel('l','sp');

s2trkern = kernel([kernel('l','s');kernel('l','sp')]);
skern_hq = kernel('hq','s',zk,kappa,d);


pxy2u = skern.eval(pxys,struct('r',rup));
pxy2d = skern.eval(pxys,struct('r',rdw));

% srcs = (-1:1).*[d;0]; alphas = exp(1i*kappa*d * (-1:1).');
srcs = (-2:2).*[d;0]; alphas = exp(1i*kappa*d * (-2:2).');
srcinfo = []; srcinfo.r = srcs; srcinfo.n = [1;0] + 0*srcs;

src2u = skern.eval(srcinfo,struct('r',rup))*alphas;
src2d = skern.eval(srcinfo,struct('r',rdw))*alphas;
% src2l = skern.eval(srcinfo,struct('r',rlf))*alphas;
% src2r = skern.eval(srcinfo,struct('r',rrt))*alphas;

src2wall = s2trkern.eval(srcinfo,rwall)*alphas;
pxy2wall = s2trkern.eval(pxys,rwall);


Gu = skern_hq.eval(struct('r',[0;0]),struct('r',rup));
Gd = skern_hq.eval(struct('r',[0;0]),struct('r',rdw));

id_l = 1:(size(rwall.r,2));
id_r = (size(rwall.r,2)+1):(2*size(rwall.r,2));

% A = [pxy2u;pxy2d;pxy2r - exp(1i*kappa*d)*pxy2l];
A = [pxy2u;pxy2d;pxy2wall(id_r,:) - exp(1i*kappa*d)*pxy2wall(id_l,:)];

% r = [Gu - src2u; Gd - src2d; -src2r + exp(1i*kappa*d)*src2l];
r = [Gu - src2u; Gd - src2d; -src2wall(id_r) + exp(1i*kappa*d)*src2wall(id_l)];

cs = A\r;


srcall = [srcs, pxy]; csall = [alphas;cs];
figure(2); clf
scatter(rup(1,:), rup(2,:),2,'.')
hold on
scatter(rdw(1,:), rdw(2,:),2,'.')
scatter(rlf(1,:), rlf(2,:),2,'.')
scatter(rrt(1,:), rrt(2,:),2,'.')
scatter(pxy(1,:), pxy(2,:),100,'.')
scatter(srcs(1,:), srcs(2,:),100,'.')
hold off
axis equal
%%
nplot = 100;
XX = linspace(-1.3*d/2,1.3*d/2,nplot);
YY = linspace(-1.3*3*d/2,1.3*3*d/2,nplot);
[X,Y] = meshgrid(XX,YY);

targs = [X(:).';Y(:).'];

% us = skern.eval(struct('r',srcall),struct('r',targs))*csall;
us = skern.eval(srcinfo,struct('r',targs))*alphas +  skern.eval(pxys,struct('r',targs))*cs;



% us_u = skern.eval(struct('r',srcall),struct('r',rup))*csall;
% us_u = skern.eval(srcinfo,struct('r',rup))*alphas +  ckern.eval(pxys,struct('r',rup))*cs;
% norm(us_u - Gu)
% 
% % us_d = skern.eval(struct('r',srcall),struct('r',rdw))*csall;
% us_u = skern.eval(srcinfo,struct('r',rdw))*alphas +  ckern.eval(pxys,struct('r',rdw))*cs;
% norm(us_d - Gd)
% 
% us_l = skern.eval(struct('r',srcall),struct('r',rlf))*csall;
% us_r = skern.eval(struct('r',srcall),struct('r',rrt))*csall;
% norm(us_r - exp(1i*kappa*d)*us_l)
% 
% 
% us_l= skern.eval(srcinfo,struct('r',rlf))*alphas + ckern.eval(pxys,struct('r',rlf))*cs;
% us_r= skern.eval(srcinfo,struct('r',rrt))*alphas + ckern.eval(pxys,struct('r',rrt))*cs;
% norm(us_r - exp(1i*kappa*d)*us_l)
% 
% 
% us_w = s2trkern.eval(srcinfo,rwall)*alphas + c2trkern.eval(pxys,rwall)*cs;
% norm(us_w(id_r) - exp(1i*kappa*d)*us_w(id_l))
% 
% us_l = skern_hq.eval(struct('r',srcall),struct('r',rlf));
% us_r = skern_hq.eval(struct('r',srcall),struct('r',rrt));
% norm(us_r - exp(1i*kappa*d)*us_l)

us_0 = skern_hq.eval(struct('r',[0;0]),struct('r',targs));

%%
figure(3);
subplot(1,3,1)
h = pcolor(X,Y,reshape(abs(us),size(X))); h.EdgeColor = 'None';
colorbar

subplot(1,3,2)
h = pcolor(X,Y,reshape(abs(us_0),size(X))); h.EdgeColor = 'None';
colorbar

subplot(1,3,3)
h = pcolor(X,Y,reshape(log10(abs(us-us_0)),size(X))); h.EdgeColor = 'None';
% h = pcolor(X,Y,reshape(abs(abs(us)-abs(us_0)),size(X))); h.EdgeColor = 'None';
hold on
scatter(rup(1,:), rup(2,:),2,'.')
scatter(rdw(1,:), rdw(2,:),2,'.')
scatter(rlf(1,:), rlf(2,:),2,'.')
scatter(rrt(1,:), rrt(2,:),2,'.')
scatter(srcall(1,:), srcall(2,:),2,'.')
hold off
colorbar


norm(us-us_0)/norm(us_0)