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
sn = chnk.helm2dquas.latticecoefs(ns,zk,d,xi,(exp(1i*xi*d)),a,M,l+1);

ising = 1;

d1 = [-1/2 0 1/2]/h;

src = []; src.r = [0.1;0.2]; 
targ = []; targ.r = [0.4;0.1] + h*[-1:1;0*(-1:1)];
[val,grad,hess,third,fourth] = chnk.flex2dquas.green(src.r,targ.r,zk,xi,d,sn,l,ising);

err1 = abs(d1*val - grad(2,1,1))  / abs(grad(2,1,1)) % G_{x}
err2 = abs(d1*grad(:,:,1) - hess(2,1,1))  / abs(hess(2,1,1)) % G_{xx}
err3 = abs(d1*hess(:,:,1) - third(2,1,1))  / abs(third(2,1,1)) % G_{xxx}
err4 = abs(d1*third(:,:,1) - fourth(2,1,1))  / abs(fourth(2,1,1)) % G_{xxxx}


targ = []; targ.r = [0.4;0.1] + h*[0*(-1:1);(-1:1)];
[val,grad,hess,third,fourth] = chnk.flex2dquas.green(src.r,targ.r,zk,xi,d,sn,l,ising);

err5 = abs(d1*val - grad(2,1,2))  / abs(grad(2,1,1)) % G_{y}
err6 = abs(d1*grad(:,:,1) - hess(2,1,2))  / abs(hess(2,1,2)) % G_{xy}
err7 = abs(d1*grad(:,:,2) - hess(2,1,3))  / abs(hess(2,1,3)) % G_{yy}

err8 = abs(d1*hess(:,:,1) - third(2,1,2))  / abs(third(2,1,2)) % G_{xxy}
err9 = abs(d1*hess(:,:,2) - third(2,1,3))  / abs(third(2,1,3)) % G_{xyy}
err10 = abs(d1*hess(:,:,3) - third(2,1,4))  / abs(third(2,1,4)) % G_{yyy}

err11 = abs(d1*third(:,:,1) - fourth(2,1,2))  / abs(fourth(2,1,2)) % G_{xxxy}
err12 = abs(d1*third(:,:,2) - fourth(2,1,3))  / abs(fourth(2,1,3)) % G_{xxyy}
err13 = abs(d1*third(:,:,3) - fourth(2,1,4))  / abs(fourth(2,1,4)) % G_{xyyy}
err14 = abs(d1*third(:,:,4) - fourth(2,1,5))  / abs(fourth(2,1,5)) % G_{yyyy}




%%
src2 = [0;0];

nplot = 100;
XX = linspace(-3*d/2,3*d/2,nplot);
% YY = linspace(-1.3*3*d/2,1.3*3*d/2,nplot);
[X,Y] = meshgrid(XX,XX);

targ2 = [X(:).';Y(:).'];
% [val2,grad2,hess2,third2,fourth2] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,ising);
[val2] = chnk.flex2dquas.green(src2,targ2,zk,xi,d,sn,l,ising);

figure(1);clf
h = pcolor(X,Y,reshape(real(val2(:,1)),size(X))); h.EdgeColor = 'None';


