d = 1;
zk = 1.5;

nnode = 80;

ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

xi = ts - 0.3i*sin(ts*d);
xi = xi(57);

h = 0.000001;

l=2; N = 80; a = 15; M = 1e4;
ns = (0:N).';
sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,xi,(exp(1i*xi*d)),a,M,l+1);
sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,xi,(exp(1i*xi*d)),a,M,l+1);
sn = cat(3,sn1,sn2);
ising = 1;

d1 = [-1/2 0 1/2]/h;

src = []; src.r = [0.1;0.2]; 
targ = []; targ.r = [0.4;0.1] + h*[-1:1;0*(-1:1)];
[val,grad,hess,third,fourth,fifth,sixth] = chnk.flex2dquas.green(src.r,targ.r,zk,xi,d,sn,l,ising);

err1 = abs(d1*val - grad(2,1,1))  / abs(grad(2,1,1)) % G_{x}
err2 = abs(d1*grad(:,:,1) - hess(2,1,1))  / abs(hess(2,1,1)) % G_{xx}
err3 = abs(d1*hess(:,:,1) - third(2,1,1))  / abs(third(2,1,1)) % G_{xxx}
err4 = abs(d1*third(:,:,1) - fourth(2,1,1))  / abs(fourth(2,1,1)) % G_{xxxx}
err5 = abs(d1*fourth(:,:,1) - fifth(2,1,1))  / abs(fifth(2,1,1)) % G_{xxxxx}
err6 = abs(d1*fifth(:,:,1) - sixth(2,1,1))  / abs(sixth(2,1,1)) % G_{xxxxxx}


targ = []; targ.r = [0.4;0.1] + h*[0*(-1:1);(-1:1)];
[val,grad,hess,third,fourth,fifth,sixth] = chnk.flex2dquas.green(src.r,targ.r,zk,xi,d,sn,l,ising);

err7 = abs(d1*val - grad(2,1,2)) % / abs(grad(2,1,1)) % G_{y}
err8 = abs(d1*grad(:,:,1) - hess(2,1,2)) % / abs(hess(2,1,2)) % G_{xy}
err9 = abs(d1*grad(:,:,2) - hess(2,1,3)) % / abs(hess(2,1,3)) % G_{yy}

err10 = abs(d1*hess(:,:,1) - third(2,1,2))  % / abs(third(2,1,2)) % G_{xxy}
err11 = abs(d1*hess(:,:,2) - third(2,1,3))  % / abs(third(2,1,3)) % G_{xyy}
err12 = abs(d1*hess(:,:,3) - third(2,1,4))  % / abs(third(2,1,4)) % G_{yyy}

err13 = abs(d1*third(:,:,1) - fourth(2,1,2))  % / abs(fourth(2,1,2)) % G_{xxxy}
err14 = abs(d1*third(:,:,2) - fourth(2,1,3))  % / abs(fourth(2,1,3)) % G_{xxyy}
err15 = abs(d1*third(:,:,3) - fourth(2,1,4))  % / abs(fourth(2,1,4)) % G_{xyyy}
err16 = abs(d1*third(:,:,4) - fourth(2,1,5))  % / abs(fourth(2,1,5)) % G_{yyyy}

err17 = abs(d1*fourth(:,:,1) - fifth(2,1,2))  % / abs(fifth(2,1,2)) % G_{xxxxy}
err18 = abs(d1*fourth(:,:,2) - fifth(2,1,3))  % / abs(fifth(2,1,3)) % G_{xxxyy}
err19 = abs(d1*fourth(:,:,3) - fifth(2,1,4))  % / abs(fifth(2,1,4)) % G_{xxyyy}
err20 = abs(d1*fourth(:,:,4) - fifth(2,1,5))  % / abs(fifth(2,1,5)) % G_{xyyyy}
err21 = abs(d1*fourth(:,:,5) - fifth(2,1,6))  % / abs(fifth(2,1,5)) % G_{yyyyy}

err22 = abs(d1*fifth(:,:,1) - sixth(2,1,2))  % / abs(sixth(2,1,2)) % G_{xxxxxy}
err23 = abs(d1*fifth(:,:,2) - sixth(2,1,3))  % / abs(sixth(2,1,3)) % G_{xxxxyy}
err24 = abs(d1*fifth(:,:,3) - sixth(2,1,4))  % / abs(sixth(2,1,4)) % G_{xxxyyy}
err25 = abs(d1*fifth(:,:,4) - sixth(2,1,5))  % / abs(sixth(2,1,5)) % G_{xxyyyy}
err26 = abs(d1*fifth(:,:,5) - sixth(2,1,6))  % / abs(sixth(2,1,6)) % G_{xyyyyy}
err27 = abs(d1*fifth(:,:,6) - sixth(2,1,7))  % / abs(sixth(2,1,7)) % G_{yyyyyy}


assert(err1 < 1e-8)
assert(err2 < 1e-8)
assert(err3 < 1e-8)
assert(err4 < 1e-8)
assert(err5 < 5e-6)
assert(err6 < 1e-4)
assert(err7 < 1e-8)
assert(err8 < 1e-8)
assert(err9 < 1e-8)
assert(err10 < 1e-8)
assert(err11 < 1e-8)
assert(err12 < 1e-8)
assert(err13 < 1e-8)
assert(err14 < 1e-8)
assert(err15 < 1e-8)
assert(err16 < 1e-8)
assert(err17 < 5e-6)
assert(err18 < 1e-6)
assert(err19 < 1e-7)
assert(err20 < 1e-7)
assert(err21 < 1e-8)
assert(err22 < 1e-3)
assert(err23 < 1e-5)
assert(err24 < 1e-5)
assert(err25 < 5e-5)
assert(err26 < 1e-6)
assert(err27 < 1e-5)



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
[val2,grad2,hess2,third2,fourth2,fifth2,sixth2] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,ising);
% [val2,grad2,hess2] = chnk.helm2dquas.green(src2,targ2,zk,xi,d,sn(:,:,1),l,ising);


[val_true, grad_true, hess_true, third_true, fourth_true] =  quasi_flex_dual_sum(X(:).',Y(:).',zk,xi,d);


%%

norm(val2(:) -val_true(:))

figure(1);clf
h = pcolor(X,Y,reshape(real(sixth2(:,:,1)),size(X))); h.EdgeColor = 'None';

%% derivatives of hkdiffgreen

zk = 1.5;
h = 0.00001;

d1 = [-1/2 0 1/2]/h;

src = []; src.r = [2.1;2.2]; 
targ = []; targ.r = [0.4;0.1] + h*[-1:1;0*(-1:1)];
[val,grad,hess,third,fourth,fifth,sixth] = chnk.flex2d.hkdiffgreen(zk,src.r,targ.r);

err1 = abs(d1*val - grad(2,1,1))  / abs(grad(2,1,1)) % G_{x}
err2 = abs(d1*grad(:,:,1) - hess(2,1,1))  / abs(hess(2,1,1)) % G_{xx}
err3 = abs(d1*hess(:,:,1) - third(2,1,1))  / abs(third(2,1,1)) % G_{xxx}
err4 = abs(d1*third(:,:,1) - fourth(2,1,1))  / abs(fourth(2,1,1)) % G_{xxxx}
err5 = abs(d1*fourth(:,:,1) - fifth(2,1,1))  / abs(fifth(2,1,1)) % G_{xxxx}
err6 = abs(d1*fifth(:,:,1) - sixth(2,1,1))  / abs(sixth(2,1,1)) % G_{xxxx}

targ = []; targ.r = [0.4;0.1] + h*[0*(-1:1);(-1:1)];
[val,grad,hess,third,fourth,fifth] = chnk.flex2d.hkdiffgreen(zk,src.r,targ.r);

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

err21 = abs(d1*fifth(:,:,1) - sixth(2,1,2))  / abs(sixth(2,1,2)) % G_{xxxxy}
err22 = abs(d1*fifth(:,:,2) - sixth(2,1,3))  / abs(sixth(2,1,3)) % G_{xxxyy}
err23 = abs(d1*fifth(:,:,3) - sixth(2,1,4))  / abs(sixth(2,1,4)) % G_{xxyyy}
err24 = abs(d1*fifth(:,:,4) - sixth(2,1,5))  / abs(sixth(2,1,5)) % G_{xyyyy}
err25 = abs(d1*fifth(:,:,5) - sixth(2,1,6))  / abs(sixth(2,1,6)) % G_{yyyyy}
err26 = abs(d1*fifth(:,:,6) - sixth(2,1,7))  / abs(sixth(2,1,7)) % G_{yyyyy}
