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

ht = 1.02*d; hb = -1.02*d;
[pxys_l, cs_l] = build_pxys(zk,xi,d,ht,hb,skern,s2trkern,l,40);

ising = 0;
fkern1 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate',xi,d,sn,pxys_l,cs_l,l,ising,nu);
fkern2 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval',xi,d,sn,pxys_l,cs_l,l,ising,nu);
hilb = @(s,t) chnk.lap2dquas.kern(s,t,'hilb',xi,d,pxys_l,cs_l,l,ising);
hilbprime = @(s,t) chnk.lap2dquas.kern(s,t,'hilbprime',xi,d,pxys_l,cs_l,l,ising);

chnkr = chunkerfuncuni(@(t) ellipse(t,2,1),32);
chnkr = chnkr.sort();

targind = 24;
kappa = signed_curvature(chnkr);

targ.r = chnkr.r(:,targind); targ.n = chnkr.n(:,targind); 
targ.d = chnkr.d(:,targind); targ.d2 = chnkr.d2(:,targind); 

src.r = chnkr.r(:,targind+10); src.n = chnkr.n(:,targind+10); 
src.d = chnkr.d(:,targind+10); src.d2 = chnkr.d2(:,targind+10); 

ref = fkern1(src,targ);
hilbref = hilb(src,targ);
hilbprimeref = hilbprime(src,targ);

k11ref = ref(1,1);
k21ref = ref(2,1) + beta/2*hilbprimeref;
k12ref = ref(1,2);
k22ref = ref(2,2);
k11href = ref(3,1) - beta.^2/2*hilbref;
k21href = ref(4,1); 

% bdry_pt = chnkr.r(:,targind);
% bdry_n = chnkr.n(:,targind);
% bdry_tau = chnkr.d(:,targind) / norm(chnkr.d(:,targind));

h = 0.001;
d2dn2 = [	15/4	-77/6	107/6	-13	61/12	-5/6] / h^2;
d2dtau2 = [-1/560	8/315	-1/5	8/5	-205/72	8/5	-1/5	8/315	-1/560] / h^2;
d3dn3 = [-17/4	71/4    -59/2	49/2	-41/4	7/4]/h^3;
ddn = [-137/60	5	-5	10/3	-5/4	1/5]/h;
ddtau = [1/280	-4/105	1/5	-4/5	0	4/5	-1/5	4/105	-1/280] / h;
d3dndtau2 = d2dtau2.*ddn.';

targ.r = chnkr.r(:,targind) +  h*(0:5).*chnkr.n(:,targind);
n_sol = fkern2(src,targ);

targ.r =  chnkr.r(:,targind) +  h*(-4:4).*chnkr.d(:,targind) / vecnorm(chnkr.d(:,targind));
tau_sol = fkern2(src,targ);
h_sol = hilb(src,targ);

k1fd = d2dn2*n_sol + nu*d2dtau2*tau_sol;

errfirstrow = k1fd - [k11ref,k11href,k12ref]

theta = atan2(chnkr.n(2,targind),chnkr.n(1,targind))-pi/2;
[xpts,ypts] = meshgrid(-4*h:h:4*h,0:h:5*h);
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];
rot_pts = R*[xpts(:) ypts(:)].';
eval_pts = rot_pts + chnkr.r(:,targind);


if (true)
figure(2);
hold on
plot(eval_pts(1,:),eval_pts(2,:),'x')
hold on
plot(chnkr,'x-')
hold on
plot(chnkr.r(1,targind),chnkr.r(2,targind))
title('Finite difference stencil')
axis square
axis equal
end

targ.r = eval_pts;
usol = fkern2(src,targ);
gnsol = reshape(usol(:,1),size(xpts));
gtsol = reshape(usol(:,2),size(xpts));
gsol = reshape(usol(:,3),size(xpts));

k21fd = d3dn3*gnsol(:,5) + (2-nu)*sum(d3dndtau2.*gnsol,'all') ...
    + (1-nu)*kappa(targind)*(d2dtau2*gnsol(1,:).' - d2dn2*gnsol(:,5)) ;

k21hfd = d3dn3*gtsol(:,5) + (2-nu)*sum(d3dndtau2.*gtsol,'all') ...
    + (1-nu)*kappa(targind)*(d2dtau2*gtsol(1,:).' - d2dn2*gtsol(:,5));

k22fd = d3dn3*gsol(:,5) + (2-nu)*sum(d3dndtau2.*gsol,'all') ...
    + (1-nu)*kappa(targind)*(d2dtau2*gsol(1,:).' - d2dn2*gsol(:,5)) ;

errsecondrow = [k21ref k21href k22ref] - [k21fd k21hfd k22fd]

errhilb = ddtau*h_sol - hilbprimeref