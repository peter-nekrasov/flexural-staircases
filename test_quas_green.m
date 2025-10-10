zk = pi/3;
d = 1;

nnode = 400;

ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

xi = ts - 0.3i*sin(ts*d);
xip = 1 - 0.3i*d*cos(ts*d);

plot(real(xi),imag(xi),'x-')

src = []; src.r = [0.1;0.2];
targ = []; targ.r = [0.4;0.9];

kern1 = kernel('hq','s',zk,xi,d);
val1 = (xip*kern1.eval(src,targ))*ws;

kern2 = kernel('h','s',zk);
val2 = kern2.eval(src,targ);

err1 = abs(val1-val2) / abs(val1);

h0qkern = kernel('hq','s',zk,xi,d);
k0qkern = kernel('hq','s',1i*zk,xi,d);

flex_kern_q = @(s,t) 1/(2*zk^2)*(h0qkern.eval(s,t) - k0qkern.eval(s,t));

val1 = xip*flex_kern_q(src,targ)*ws;

flex_kern = @(s,t) chnk.flex2d.kern(zk,s,t,'s');

val2 = flex_kern(src,targ);

abs(val1-val2)