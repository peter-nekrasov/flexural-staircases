d = 1;

kappa = 0.001;

ht = 1.02*d; hb = -1.02*d;
l = 2;
npxy = 40;

skern = kernel('l','s');
s2trkern = kernel([kernel('l','s');kernel('l','sp')]);


tic;
[pxys, cs] = build_pxys(0,kappa,d,ht,hb,skern,s2trkern,l,npxy);
toc;


src = []; src.r = [0;3]; src.n = [1;0];

nplot = 100;
XX = linspace(-3*d/2,3*d/2,nplot);
[X,Y] = meshgrid(XX,XX);

targ = []; 
targ.r = [X(:).';Y(:).'];
targ.n = [1;1]/sqrt(2) + 0*targ.r;


u = quasi_lap_kern(src,targ,'d',kappa,d,pxys,cs,l);


 figure(1);clf
h = pcolor(X,Y,reshape(real(u),size(X))); h.EdgeColor = 'None';
colorbar


chnkr = chunkerfunc(@starfish,struct('eps',1e-10)); chnkr = 0.25*chnkr;
chnkrs = [];
for i = -5:5
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

figure(2);clf
plot(chnkr)
hold on
plot(chnkrs,'k')
hold off
axis equal


dkern = kernel(@(s,t) quasi_lap_kern(s,t,'d',kappa,d,pxys,cs,l));
dkern.sing = 'smooth';
hkern = kernel(@(s,t) quasi_lap_kern(s,t,'hilb',kappa,d,pxys,cs,l));
hkern.sing = 'pv';

dmat = chunkermat(chnkr,dkern);
hmat = chunkermat(chnkr,hkern);

rhs = quasi_lap_kern(src,chnkr,'s',kappa,d,pxys,cs,l);
a = 0.25*hmat * (hmat * rhs);
b = -0.25*rhs + dmat * (dmat * rhs);

norm(a-b) / norm(a)



dkern = kernel(@(s,t) chnk.lap2d.kern(s,t,'d'));
dkern.sing = 'smooth';
hkern = kernel(@(s,t) chnk.lap2d.kern(s,t,'hilb'));
hkern.sing = 'pv';

dmat_0 = chunkermat(chnkr,dkern);
hmat_0 = chunkermat(chnkr,hkern);


% norm(0.25 * hmat_0*hmat_0 - (-0.25*eye(size(hmat_0)) + dmat_0*dmat_0))

rhs = skern.eval(src,chnkr);

a = 0.25*hmat_0 * (hmat_0 * rhs);
b = -0.25*rhs + dmat_0 * (dmat_0 * rhs);

norm(a-b) / norm(a)
%%

cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
nch = 20; A = 1;
chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
chnkr = reverse(chnkr);

chnkrs = [];
for i = -5:5
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

figure(2);clf
plot(chnkr,'.')
hold on
plot(chnkrs,'k.')
scatter(src.r(1,:),src.r(2,:))
hold off
axis equal


dkern = kernel(@(s,t) quasi_lap_kern(s,t,'d',kappa,d,pxys,cs,l));
dkern.sing = 'smooth';
hkern = kernel(@(s,t) quasi_lap_kern(s,t,'hilb',kappa,d,pxys,cs,l));
hkern.sing = 'pv';

dmat = chunkermat(chnkr,dkern);
hmat = chunkermat(chnkr,hkern);

rhs = quasi_lap_kern(src,chnkr,'s',kappa,d,pxys,cs,l);

a = 0.25*hmat * (hmat * rhs);
b = -0.25*rhs + dmat * (dmat * rhs);

norm(a-b) / norm(a)

function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end