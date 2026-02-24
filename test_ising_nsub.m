%%
d = 1.2;
src2 = [0.2;0.2];

zk = 4;
xi = 0.5;

alpha = exp(1i*xi*d);

nplot = 100;
XX = linspace(-5*d/2,5*d/2,nplot);
% YY = linspace(-1.3*3*d/2,1.3*3*d/2,nplot);
[X,Y] = meshgrid(XX,XX);

targ2 = [X(:).';Y(:).'];

l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,xi,(exp(1i*xi*d)),a,M,l+1);
[s0_l,sn_l] = chnk.lap2dquas.latticecoefs((1:N),d,xi,l);
%

[val1] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,1);


[val2] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,0,0);
val2 = val2 + chnk.flex2d.hkdiffgreen(zk,src2,targ2);

[val2] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,0,1);
val2 = val2 + chnk.flex2d.hkdiffgreen(zk,src2,targ2) ...
    + alpha.^(-1).*chnk.flex2d.hkdiffgreen(zk,src2 + [-d;0],targ2) ...
    + alpha.*chnk.flex2d.hkdiffgreen(zk,src2 + [d;0],targ2);

figure(1);clf
h = pcolor(X,Y,reshape(log10(abs(val1 - val2) ./ max(abs(val1(:)))),size(X))); h.EdgeColor = 'None'; colorbar
