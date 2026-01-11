d = 1;
zk = 1.5;

nnode = 80;

ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

xi = ts - 0.3i*sin(ts*d);
xi = xi(57);

h = 0.00001;

l=2; N = 40; a = 15; M = 1e4;
ns = (0:N).';
sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,xi,(exp(1i*xi*d)),a,M,l+1);
sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,xi,(exp(1i*xi*d)),a,M,l+1);
sn = cat(3,sn1,sn2);
ising = 1;

d1 = [-1/2 0 1/2]/h;

src = []; src.r = [0.1;0.2]; 
targ = []; targ.r = [0.4;0.1] + h*[-1:1;0*(-1:1)];
[val,grad,hess,third,fourth,fifth] = chnk.flex2dquas.green(src.r,targ.r,zk,xi,d,sn,l,ising);

err1 = abs(d1*val - grad(2,1,1))  / abs(grad(2,1,1)) % G_{x}
err2 = abs(d1*grad(:,:,1) - hess(2,1,1))  / abs(hess(2,1,1)) % G_{xx}
err3 = abs(d1*hess(:,:,1) - third(2,1,1))  / abs(third(2,1,1)) % G_{xxx}
err4 = abs(d1*third(:,:,1) - fourth(2,1,1))  / abs(fourth(2,1,1)) % G_{xxxx}
err5 = abs(d1*fourth(:,:,1) - fifth(2,1,1))  / abs(fifth(2,1,1)) % G_{xxxx}

targ = []; targ.r = [0.4;0.1] + h*[0*(-1:1);(-1:1)];
[val,grad,hess,third,fourth,fifth] = chnk.flex2dquas.green(src.r,targ.r,zk,xi,d,sn,l,ising);

err6 = abs(d1*val - grad(2,1,2))  / abs(grad(2,1,1)) % G_{y}
err7 = abs(d1*grad(:,:,1) - hess(2,1,2))  / abs(hess(2,1,2)) % G_{xy}
err8 = abs(d1*grad(:,:,2) - hess(2,1,3))  / abs(hess(2,1,3)) % G_{yy}

err9 = abs(d1*hess(:,:,1) - third(2,1,2))  / abs(third(2,1,2)) % G_{xxy}
err10 = abs(d1*hess(:,:,2) - third(2,1,3))  / abs(third(2,1,3)) % G_{xyy}
err11 = abs(d1*hess(:,:,3) - third(2,1,4))  / abs(third(2,1,4)) % G_{yyy}

err12 = abs(d1*third(:,:,1) - fourth(2,1,2))  / abs(fourth(2,1,2)) % G_{xxxy}
err13 = abs(d1*third(:,:,2) - fourth(2,1,3))  / abs(fourth(2,1,3)) % G_{xxyy}
err14 = abs(d1*third(:,:,3) - fourth(2,1,4))  / abs(fourth(2,1,4)) % G_{xyyy}
err15 = abs(d1*third(:,:,4) - fourth(2,1,5))  / abs(fourth(2,1,5)) % G_{yyyy}

err16 = abs(d1*fourth(:,:,1) - fifth(2,1,2))  / abs(fifth(2,1,2)) % G_{xxxxy}
err17 = abs(d1*fourth(:,:,2) - fifth(2,1,3))  / abs(fifth(2,1,3)) % G_{xxxyy}
err18 = abs(d1*fourth(:,:,3) - fifth(2,1,4))  / abs(fifth(2,1,4)) % G_{xxyyy}
err19 = abs(d1*fourth(:,:,4) - fifth(2,1,5))  / abs(fifth(2,1,5)) % G_{xyyyy}
err20 = abs(d1*fourth(:,:,5) - fifth(2,1,6))  / abs(fifth(2,1,5)) % G_{yyyyy}


assert(err1 < 1e-8)
assert(err2 < 1e-8)
assert(err3 < 1e-8)
assert(err4 < 1e-8)
assert(err5 < 5e-6)
assert(err6 < 1e-8)
assert(err7 < 1e-8)
assert(err8 < 1e-8)
assert(err9 < 1e-8)
assert(err10 < 1e-8)
assert(err11 < 1e-8)
assert(err12 < 1e-8)
assert(err13 < 1e-8)
assert(err14 < 1e-8)
assert(err15 < 1e-8)
assert(err16 < 5e-6)
assert(err17 < 1e-6)
assert(err18 < 1e-7)
assert(err19 < 1e-7)
assert(err20 < 1e-7)


return

%%
src2 = [0;0];

nplot = 100;
XX = linspace(-3*d/2,3*d/2,nplot);
% YY = linspace(-1.3*3*d/2,1.3*3*d/2,nplot);
[X,Y] = meshgrid(XX,XX);

targ2 = [X(:).';Y(:).'];

% xi1 = xi(57);
% sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,xi1,(exp(1i*xi*d)),a,M,l+1);
% sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,xi1,(exp(1i*xi*d)),a,M,l+1);
% sn = cat(3,sn1,sn2);

% [val2,grad2,hess2,third2,fourth2] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,ising);
% [val2] = chnk.helm2dquas.green(src2,targ2,zk,xi,d,sn(:,:,1),l,ising);
[val2,grad2,hess2,third2,fourth2] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,ising);
% [val2,grad2,hess2] = chnk.helm2dquas.green(src2,targ2,zk,xi,d,sn(:,:,1),l,ising);


[val_true, grad_true, hess_true, third_true, fourth_true] =  quasi_flex_dual_sum(X(:).',Y(:).',zk,xi,d);




norm(val2(:) -val_true(:))

figure(1);clf
h = pcolor(X,Y,reshape(real(val2(:,:,1)),size(X))); h.EdgeColor = 'None';

