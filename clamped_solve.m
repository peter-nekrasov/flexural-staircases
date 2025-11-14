
zk = 1.2;
d = 1.2;

kappa = 0.5;

src = []; src.r = [0;-2]; src.n = [1;0];

nplot = 80;
xx = linspace(-1.5*d, 1.5*d,nplot);
yy = xx;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];


if false
    cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
    nch = 20; A = 1;
    chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
    chnkr = reverse(chnkr);
    wtarg = cos_func(targ.r(1,:),d,A) ;
    iout = targ.r(2,:) > wtarg(2,:);
else
    chnkr = chunkerfunc(@starfish,struct('eps',1e-10)); chnkr = 0.25*chnkr;
    targmod = []; targmod.r = [targ.r(1,:), targ.r(2,:)];
    iout = ~chunkerinterior(chnkr,targmod);
end

targout = []; targout.r = targ.r(:,iout);

chnkrs = [];
for i = -5:5
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = 1;

figure(2);clf
quiver(chnkr)
hold on
plot(chnkrs,'k.')
scatter(src.r(1,:),src.r(2,:))
h = pcolor(X,Y, reshape(real(us),size(X))); h.EdgeColor = 'None';
hold off
axis equal





function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end