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

[val_ref,grad_ref,hess_ref,third_ref,fourth_ref,fifth_ref,sixth_ref] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,1);

nsub = 2;
[val,grad,hess,third,fourth,fifth,sixth] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,0,nsub);

for ii = -nsub:nsub

[val0,grad0,hess0,third0,fourth0,fifth0,sixth0] = chnk.flex2d.hkdiffgreen(zk,src2+ii*[d;0],targ2);

val = val + alpha.^ii.*val0;
grad = grad + alpha.^ii.*grad0;
hess = hess + alpha.^ii.*hess0;
third = third + alpha.^ii.*third0;
fourth = fourth + alpha.^ii.*fourth0;
fifth = fifth + alpha.^ii.*fifth0;
sixth = sixth + alpha.^ii.*sixth0;

end

figure(1);clf
h = pcolor(X,Y,reshape(log10(abs(sixth_ref(:,:,1) - sixth(:,:,1)) ./ max(abs(fifth(:)))),size(X))); h.EdgeColor = 'None'; colorbar
