d = 1;
zk = 1.5;
nu = 0.3;
beta = (1+nu)/2;

nnode = 80;

ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

xi = 1.4-0.3i;

l=2; N = 40; a = 15; M = 1e4;
ns = (0:N).';
sn1 = chnk.helm2dquas.latticecoefs(ns,zk,d,xi,(exp(1i*xi*d)),a,M,l+1);
sn2 = chnk.helm2dquas.latticecoefs(ns,1i*zk,d,xi,(exp(1i*xi*d)),a,M,l+1);
sn = cat(3,sn1,sn2);

skern = kernel('l','s');
s2trkern = kernel([kernel('l','s');kernel('l','sp')]);

ising = 0;
fkern1 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate',xi,d,sn,[],[],l,ising,nu);
fkern2 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_eval',xi,d,sn,[],[],l,ising,nu);

chnkr = chunkerfuncuni(@(t) ellipse(t,2,1),32);
chnkr = chnkr.sort();

targind = 24;

targ.r = chnkr.r(:,targind); targ.n = chnkr.n(:,targind); 
targ.d = chnkr.d(:,targind); targ.d2 = chnkr.d2(:,targind); 

src.r = chnkr.r(:,targind+10); src.n = chnkr.n(:,targind+10); 
src.d = chnkr.d(:,targind+10); src.d2 = chnkr.d2(:,targind+10); 
curv = signed_curvature(chnkr);
kp = arclengthder(chnkr,curv);
kpp = arclengthder(chnkr,kp);

src.data = zeros(2,1);
src.data(1) = kp(targind+10);
src.data(2) = kpp(targind+10);

ref = fkern1(src,targ);

k11ref = ref(1,1);
k21ref = ref(2,1);
k12ref = ref(1,2);
k22ref = ref(2,2);

h = 0.001;
d2dn2 = [	15/4	-77/6	107/6	-13	61/12	-5/6] / h^2;
d2dtau2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560] / h^2;

targ.r = chnkr.r(:,targind) +  h*(0:5).*chnkr.n(:,targind);
n_sol = fkern2(src,targ);

targ.r =  chnkr.r(:,targind) +  h*(-4:4).*chnkr.d(:,targind) / vecnorm(chnkr.d(:,targind));
tau_sol = fkern2(src,targ);

targ.r =  chnkr.r(:,targind);
k1fd = fkern2(src,targ);

k2fd = d2dn2*n_sol + nu*d2dtau2*tau_sol;

err = abs(ref - [k1fd;k2fd])

%% checking bcs

fkern3 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_bcs',xi,d,sn,[],[],l,ising,nu);
fkern4 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',xi,d,sn,[],[],l,ising,nu);

ref = fkern3(src,targ);

targ.r = chnkr.r(:,targind) +  h*(0:5).*chnkr.n(:,targind);
n_sol = fkern4(src,targ);

targ.r =  chnkr.r(:,targind) +  h*(-4:4).*chnkr.d(:,targind) / vecnorm(chnkr.d(:,targind));
tau_sol = fkern4(src,targ);

targ.r =  chnkr.r(:,targind);
k1fd = fkern4(src,targ);

k2fd = d2dn2*n_sol + nu*d2dtau2*tau_sol;

err2 = abs(ref - [k1fd;k2fd])