function [pxys, cs] = build_flex_pxys(zk,kappas,d,ht,hb,l,npxy)
kappas = kappas(:).';

n = 1000;
thetas = 2*pi*(1:npxy)/npxy;
pxy0 = [cos(thetas);sin(thetas)];
% pxy = 5*d/2 * ;

pxy = [3/2 * pxy0, 2*pxy0, 3*pxy0]*d;
npxy = size(pxy,2);
pxys = []; pxys.r = pxy; pxys.n = [cos(thetas);sin(thetas)];

xs = linspace(-d/2,d/2,n);
ys = linspace(hb,ht,2*n);


rup = [xs(:).'; 0*xs(:).' + ht];
rdw = [xs(:).'; 0*xs(:).' + hb];

rwall = [];
rwall.r = [[0*ys(:).' - d/2; ys(:).'],[0*ys(:).' + d/2; ys(:).']];
% rwall.n = [1;0] + 0*rwall.r;

[pxy2u, pxy2u_grad, pxy2u_hess, pxy2u_third] = chnk.flex2d.hkdiffgreen(zk,pxys.r,rup);
dnpxy2u = pxy2u_grad(:,:,end);
dn2pxy2u = pxy2u_hess(:,:,end);
dn3pxy2u = pxy2u_third(:,:,end);

[pxy2d, pxy2d_grad, pxy2d_hess, pxy2d_third] = chnk.flex2d.hkdiffgreen(zk,pxys.r,rdw);
dnpxy2d = pxy2d_grad(:,:,end);
dn2pxy2d = pxy2d_hess(:,:,end);
dn3pxy2d = pxy2d_third(:,:,end);

[pxy2wall, pxy2wall_grad, pxy2wall_hess, pxy2wall_third] = chnk.flex2d.hkdiffgreen(zk,pxys.r,rwall.r);
dnpxy2wall = pxy2wall_grad(:,:,1);
dn2pxy2wall = pxy2wall_hess(:,:,1);
dn3pxy2wall = pxy2wall_third(:,:,1);

srcs = (-l:l).*[d;0]; alphas = exp(1i*kappas*d .* (-l:l).');
srcinfo = []; srcinfo.r = srcs; srcinfo.n = [0;0] + 0*srcs;


[src2u, src2u_grad, src2u_hess, src2u_third] = chnk.flex2d.hkdiffgreen(zk,srcs,rup);
src2u = src2u*alphas;
dnsrc2u = src2u_grad(:,:,end)*alphas;
dn2src2u = src2u_hess(:,:,end)*alphas;
dn3src2u = src2u_third(:,:,end)*alphas;

[src2d, src2d_grad, src2d_hess, src2d_third] = chnk.flex2d.hkdiffgreen(zk,srcs,rdw);
src2d = src2d*alphas;
dnsrc2d = src2d_grad(:,:,end)*alphas;
dn2src2d = src2d_hess(:,:,end)*alphas;
dn3src2d = src2d_third(:,:,end)*alphas;

[src2wall, src2wall_grad, src2wall_hess, src2wall_third] = chnk.flex2d.hkdiffgreen(zk,srcs,rwall.r);
src2wall = src2wall*alphas;
dnsrc2wall = src2wall_grad(:,:,1)*alphas;
dn2src2wall = src2wall_hess(:,:,1)*alphas;
dn3src2wall = src2wall_third(:,:,1)*alphas;


% src2u = skern.eval(srcinfo,struct('r',rup))*alphas;
% src2d = skern.eval(srcinfo,struct('r',rdw))*alphas;
% 
% src2wall = s2trkern.eval(srcinfo,rwall)*alphas;
% pxy2wall = s2trkern.eval(pxys,rwall);
%
% Gu = quasi_dual_sum(rup(1,:),rup(2,:),zk,kappas,d);
% Gd = quasi_dual_sum(rdw(1,:),rdw(2,:),zk,kappas,d);


[Gu, Gu_grad, Gu_hess, Gu_third] = quasi_flex_dual_sum(rup(1,:),rup(2,:),zk,kappas,d);
dnGu = Gu_grad(:,:,end);
dn2Gu = Gu_hess(:,:,end);
dn3Gu = Gu_third(:,:,end);

[Gd, Gd_grad, Gd_hess, Gd_third] = quasi_flex_dual_sum(rdw(1,:),rdw(2,:),zk,kappas,d);
dnGd = Gd_grad(:,:,end);
dn2Gd = Gd_hess(:,:,end);
dn3Gd = Gd_third(:,:,end);

% [Gwall, Gwall_grad, Gwall_hess, Gwall_third] = quasi_flex_dual_sum(rwall.r(1,:),rwall.r(2,:),zk,kappas,d);
% dnGwall = Gwall_grad(:,:,1);
% dn2Gwall = Gwall_hess(:,:,1);
% dn3Gwall = Gwall_third(:,:,1);

id_l = 1:(size(rwall.r,2))/2;
id_r = (size(rwall.r,2)/2+1):size(rwall.r,2);

cs = zeros(npxy, length(kappas));
for i = 1:length(kappas)
    % A = [pxy2u;pxy2d;pxy2wall(id_r,:) - exp(1i*kappas(i)*d).*pxy2wall(id_l,:)];

    A  = [pxy2u;dnpxy2u;dn2pxy2u;dn3pxy2u; ... 
        pxy2d;dnpxy2d;dn2pxy2d;dn3pxy2d; ...
        pxy2wall(id_r,:) - exp(1i*kappas(i)*d).*pxy2wall(id_l,:);
        dnpxy2wall(id_r,:) - exp(1i*kappas(i)*d).*dnpxy2wall(id_l,:);
        dn2pxy2wall(id_r,:) - exp(1i*kappas(i)*d).*dn2pxy2wall(id_l,:);
        dn3pxy2wall(id_r,:) - exp(1i*kappas(i)*d).*dn3pxy2wall(id_l,:)];

    r = [Gu(i,:).' - src2u(:,i); 
        dnGu(i,:).' - dnsrc2u(:,i); 
        dn2Gu(i,:).' - dn2src2u(:,i); 
        dn3Gu(i,:).' - dn3src2u(:,i);
        Gd(i,:).' - src2d(:,i); 
        dnGd(i,:).' - dnsrc2d(:,i); 
        dn2Gd(i,:).' - dn2src2d(:,i); 
        dn3Gd(i,:).' - dn3src2d(:,i);
        -src2wall(id_r,i) + exp(1i*kappas(i)*d).*src2wall(id_l,i);
        -dnsrc2wall(id_r,i) + exp(1i*kappas(i)*d).*dnsrc2wall(id_l,i);
        -dn2src2wall(id_r,i) + exp(1i*kappas(i)*d).*dn2src2wall(id_l,i);
        -dn3src2wall(id_r,i) + exp(1i*kappas(i)*d).*dn3src2wall(id_l,i)];

    % A  = [pxy2u;dnpxy2u;dn2pxy2u;dn3pxy2u; ... 
    %     pxy2d;dnpxy2d;dn2pxy2d;dn3pxy2d];
    % 
    % r = [Gu(i,:).' - src2u(:,i); 
    %     dnGu(i,:).' - dnsrc2u(:,i); 
    %     dn2Gu(i,:).' - dn2src2u(:,i); 
    %     dn3Gu(i,:).' - dn3src2u(:,i);
    %     Gd(i,:).' - src2d(:,i); 
    %     dnGd(i,:).' - dnsrc2d(:,i); 
    %     dn2Gd(i,:).' - dn2src2d(:,i); 
    %     dn3Gd(i,:).' - dn3src2d(:,i)];

    % A  = [pxy2u;dnpxy2u; ... 
    %     pxy2d;dnpxy2d; ...
    %     pxy2wall(id_r,:) - exp(1i*kappas(i)*d).*pxy2wall(id_l,:);
    %     dnpxy2wall(id_r,:) - exp(1i*kappas(i)*d).*dnpxy2wall(id_l,:);
    %     dn2pxy2wall(id_r,:) - exp(1i*kappas(i)*d).*dn2pxy2wall(id_l,:);
    %     dn3pxy2wall(id_r,:) - exp(1i*kappas(i)*d).*dn3pxy2wall(id_l,:)];
    % 
    % r = [Gu(i,:).' - src2u(:,i); 
    %     dnGu(i,:).' - dnsrc2u(:,i); 
    %     Gd(i,:).' - src2d(:,i); 
    %     dnGd(i,:).' - dnsrc2d(:,i); 
    %     -src2wall(id_r,i) + exp(1i*kappas(i)*d).*src2wall(id_l,i);
    %     -dnsrc2wall(id_r,i) + exp(1i*kappas(i)*d).*dnsrc2wall(id_l,i);
    %     -dn2src2wall(id_r,i) + exp(1i*kappas(i)*d).*dn2src2wall(id_l,i);
    %     -dn3src2wall(id_r,i) + exp(1i*kappas(i)*d).*dn3src2wall(id_l,i)];
    cs(:,i) = A\r;
end
pxys = pxys.r;
end