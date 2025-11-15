
zk = 1.2;
d = 1.2;

kappa = pi/d;



nplot = 80;
xx = linspace(-1.5*d, 1.5*d,nplot);
yy = xx;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];


if true
    cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
    nch = 20; A = 1;
    chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
    chnkr = reverse(chnkr);
    wtarg = cos_func(targ.r(1,:),d,A) ;
    iout = targ.r(2,:) > wtarg(2,:);
    src = []; src.r = [0;-2]; src.n = [1;0];
    % src = []; src.r = [0;2]; src.n = [1;0];
else
    chnkr = chunkerfunc(@starfish,struct('eps',1e-10)); chnkr = 0.25*chnkr;
    targmod = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
    iout = ~chunkerinterior(chnkr,targmod);
    src = []; src.r = [0;0]; src.n = [1;0];
end

targout = []; targout.r = targ.r(:,iout);

%%


skern = kernel('l','s');
s2trkern = kernel([kernel('l','s');kernel('l','sp')]);

ht = 1.02*d; hb = -1.02*d;
l = 2;

[pxys_f, cs_f] = build_flex_pxys(zk,kappa,d,ht,hb,l,100);
[pxys_l, cs_l] = build_pxys(zk,kappa,d,ht,hb,skern,s2trkern,l,40);

%%
fkern =  @(s,t) qflex_kern(zk, s, t, 'clamped_plate',kappa,d,pxys_f,cs_f,pxys_l,cs_l,l);

curv = signed_curvature(chnkr);
curv = curv(:);

opts = [];
opts.sing = 'log';

start = tic;
sys = chunkermat(chnkr,fkern, opts);
sys = sys - 0.5*eye(2*chnkr.npt);
sys(2:2:end,1:2:end) = sys(2:2:end,1:2:end) + curv.*eye(chnkr.npt);
toc(start)
%%

skern =  @(s,t) qflex_kern(zk, s, t, 's',kappa,d,pxys_f,cs_f,pxys_l,cs_l,l);
bskern =  @(s,t) qflex_kern(zk, s, t, 'clamped_plate_bcs',kappa,d,pxys_f,cs_f,pxys_l,cs_l,l);

rhs = -bskern(src,chnkr);

% Solving linear system
sol = sys\rhs;

ikern = @(s,t) qflex_kern(zk, s, t, 'clamped_plate_eval',kappa,d,pxys_f,cs_f,pxys_l,cs_l,l);

opts = []; opts.forcesmooth = true;
start1 = tic;
uscat = chunkerkerneval(chnkr, ikern,sol, targout,opts);
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)

uin = skern(src,targout);
utot = uscat(:)+uin(:);



%%
chnkrs = [];
for i = -1:1
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = utot;

figure(2);clf
quiver(chnkr)
hold on
plot(chnkrs,'k.')
scatter(src.r(1,:),src.r(2,:))
h = pcolor(X,Y, reshape(log10(abs(us)),size(X))); h.EdgeColor = 'None';
colorbar
hold off
axis equal





function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end