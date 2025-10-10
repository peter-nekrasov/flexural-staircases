% problem parameters

kappa = [pi-0.1-1i,0.2,-0.1,0.25i];
nkappa = length(kappa);

src = []; src.r = [[0;-1.1],[1;-1],[0.1;-0.3]]; 
targ = []; targ.r = [[1.1;0.3],[2;0]]; 
ns = size(src.r,2);
nt = size(targ.r,2);

h0qkern = kernel('hq','s',zk,kappa,d);
k0qkern = kernel('hq','s',1i*zk,kappa,d);

flex_kern_q = @(s,t) 1/(2*zk^2)*(h0qkern.eval(s,t) - k0qkern.eval(s,t));

sval = flex_kern_q(src,targ);

h0qkern = kernel('hq','s',zk,kappa,d,[],[],0);
k0qkern = kernel('hq','s',1i*zk,kappa,d,[],[],0);

flex_kern_q1a = @(s,t) 1/(2*zk^2)*(h0qkern.eval(s,t) - k0qkern.eval(s,t));

sval1b = flex_kern(src,targ); 
sval1b = repmat(reshape(sval1b,1,nt,ns), nkappa,1,1);
sval1 = flex_kern_q1a(src,targ) + reshape(sval1b,nt*nkappa,ns);

norm(sval-sval1)