%%
d = 1.2;

src = chnkr; targ = chnkr; % + [d;0];

zk = 4;
xi = 0.5;

alpha = exp(1i*xi*d);

nplot = 100;
XX = linspace(-5*d/2,5*d/2,nplot);
% YY = linspace(-1.3*3*d/2,1.3*3*d/2,nplot);
[X,Y] = meshgrid(XX,XX);

targ2 = []; targ2.r = [X(:).';Y(:).'];

l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,xi,(exp(1i*xi*d)),a,M,l+1);
[s0_l,sn_l] = chnk.lap2dquas.latticecoefs((1:N),d,xi,l);
%

fmat = chnk.flex2dquas.kern(zk, src, targ, 'free_plate',xi,d,sn,s0_l,sn_l,l,1,nu);
% fmat = fmat + chnk.flex2d.kern(zk, src, targ, 'free_plate',nu);

nsub = 1;
fmat1 = chnk.flex2dquas.kern(zk, src, targ, 'free_plate',xi,d,sn,s0_l,sn_l,l,0,nu,nsub);

for ii = -nsub:nsub
fmat0 = chnk.flex2d.kern(zk, src + ii*[d;0], targ, 'free_plate',nu);
fmat1 = fmat1 + alpha.^(ii).*fmat0;
end

for ii = 1:chnkr.npt
    id2 = (1:2) + (ii-1)*2;
    id1 = (1:4) + (ii-1)*4;
    fmat(id1,id2) = 0;
    fmat1(id1,id2) = 0;
end

errmat = fmat - fmat1;
vecnorm(errmat(:)) ./ vecnorm(fmat(:))

return
