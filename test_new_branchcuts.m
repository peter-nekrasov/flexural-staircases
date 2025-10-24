d = 1;
zk = 1.5;
zk = 4;


ht = 1.02*d; hb = -1.02*d;
l = 2;
npxy = 40;

if abs(zk) > 1e-10
    skern = kernel('h','s',zk);
    s2trkern = kernel([kernel('h','s',zk);kernel('h','sp',zk)]);
else
    skern = kernel('l','s');
    s2trkern = kernel([kernel('l','s');kernel('l','sp')]);
end


nnode = 80;

ts = linspace(-pi/d,pi/d,nnode)-1e-8;
ts = ts(2:end);
ws = 1/(nnode-1);


amp = (-1)^(mod(ceil((real(zk)) / (2*pi/d)+.5) ,2)) * 0.3;

xi = ts + amp*1i*sin(ts*d);
xip = 1 + amp*1i*d*cos(ts*d);

ws = ws*xip;

plot(real(xi),imag(xi),'x-')



%%
src = []; src.r = [0.1;0.2];
targ = []; targ.r = [0.4;0.1];

nplot = 80;
xx = linspace(-pi/d,pi/d,nplot);
[X,Y] = meshgrid(xx,xx);
xi_grid = X(:) + 1i*Y(:);

tic;
[pxys, cs] = build_pxys(zk,xi_grid,d,ht,hb,skern,s2trkern,l,npxy);
toc;

tic;
[val, grad, hess] = new_green(src.r,targ.r,zk,xi_grid,d,pxys,cs,l,1);
toc;
val = reshape(val, length(xi_grid), []);

%%

figure(1);clf
h = pcolor(X,Y, reshape(angle(val),size(X))); h.EdgeColor = 'None';
hold on, 
plot(xi,'r-')
hold off
colormap('hsv')
colorbar


%%


tic;
[pxys, cs] = build_pxys(zk,xi,d,ht,hb,skern,s2trkern,l,npxy);
toc;

tic;
[val, grad, hess] = new_green(src.r,targ.r,zk,xi,d,pxys,cs,l,1);
toc;

val1 = ws*reshape(val, length(xi), []);


kern2 = kernel('h','s',zk);
val2 = kern2.eval(src,targ);

err1 = abs(val1-val2) / abs(val1)
[val1,val2]


function fkern = flex_kernel(zk,xi,d)

fkern = 1/(2*zk^2)*(kernel('hq','s',zk,xi,d) - ...
    kernel('hq','s',1i*zk,xi,d));

end
