function [pxys, cs] = build_pxys(zk,kappas,d,ht,hb,skern,s2trkern,l,npxy)
kappas = kappas(:).';

n = 1000;
thetas = 2*pi*(1:npxy)/npxy;
pxy = 5*d/2 * [cos(thetas);sin(thetas)];
pxys = []; pxys.r = pxy; pxys.n = [cos(thetas);sin(thetas)];

xs = linspace(-d/2,d/2,n);
ys = linspace(hb,ht,2*n);


rup = [xs(:).'; 0*xs(:).' + ht];
rdw = [xs(:).'; 0*xs(:).' + hb];

rwall = [];
rwall.r = [[0*ys(:).' - d/2; ys(:).'],[0*ys(:).' + d/2; ys(:).']];
rwall.n = [1;0] + 0*rwall.r;

pxy2u = skern.eval(pxys,struct('r',rup));
pxy2d = skern.eval(pxys,struct('r',rdw));


srcs = (-l:l).*[d;0]; alphas = exp(1i*kappas*d .* (-l:l).');
srcinfo = []; srcinfo.r = srcs; srcinfo.n = [0;0] + 0*srcs;

src2u = skern.eval(srcinfo,struct('r',rup))*alphas;
src2d = skern.eval(srcinfo,struct('r',rdw))*alphas;

src2wall = s2trkern.eval(srcinfo,rwall)*alphas;
pxy2wall = s2trkern.eval(pxys,rwall);


Gu = quasi_dual_sum(rup(1,:),rup(2,:),zk,kappas,d);
Gd = quasi_dual_sum(rdw(1,:),rdw(2,:),zk,kappas,d);

id_l = 1:(size(rwall.r,2));
id_r = (size(rwall.r,2)+1):(2*size(rwall.r,2));

cs = zeros(npxy, length(kappas));
for i = 1:length(kappas)
    A = [pxy2u;pxy2d;pxy2wall(id_r,:) - exp(1i*kappas(i)*d).*pxy2wall(id_l,:)];
    r = [Gu(i,:).' - src2u(:,i); Gd(i,:).' - src2d(:,i); -src2wall(id_r,i) + exp(1i*kappas(i)*d).*src2wall(id_l,i)];
    cs(:,i) = A\r;
end
pxys = pxys.r;
end