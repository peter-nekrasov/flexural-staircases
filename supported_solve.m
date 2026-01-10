zk = 2;
zk = 0.2;
d = 1.2;
% d = 2*pi;
nu = 0.3;

% kappa = pi/d;
kappa = 0.5;
% kappa = pi/d;
ws = 1;

nnode = 62;
ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa = ts + amp*1i*sin(ts*d);
xip = 1 + amp*1i*d*cos(ts*d);
ws = ws*xip;

nkappa = length(kappa);

nplot = 80;
xx = linspace(-1.5*d, 1.5*d,nplot);
yy = xx;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];


if false
    cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
    cparams = []; cparams.ta = 0; cparams.tb = d;
    nch = 20/2; A = 0.2;
    chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
    chnkr = reverse(chnkr);
    wtarg = cos_func(targ.r(1,:),d,A) ;
    iout = targ.r(2,:) > wtarg(2,:);
    src = []; src.r = [[-0.2;-2],[0.2;-2]]; src.n = [[1;0],[1;0]];
    % src = []; src.r = [0;2]; src.n = [1;0];
    chnkr = makedatarows(chnkr,2);

else
    chnkr = chunkerfunc(@(t) starfish(t,3,0.1),struct('eps',1e-10,'maxchunklen',0.4)); chnkr = 0.25*chnkr;
    chnkr = chnkr*1.3;
    targmod = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
    iout = ~chunkerinterior(chnkr,targmod);
    src = []; src.r = [[-0.01;-0.01], [0.05;0.02]]; %src.n = [1;0];
    chnkr = makedatarows(chnkr,2);
end

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);
targout = []; targout.r = targmod.r(:,iout);
targout_0 = []; targout_0.r = targ.r(:,iout);


l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);

%%
ising = 0;
fkern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate',kappa,d,sn,[],[],l,ising,nu);

curv = signed_curvature(chnkr);
kp = arclengthder(chnkr,curv);
kpp = arclengthder(chnkr,kp);

% supported plate kernels expect (d/ds) kappa in the first data row
% and (d^2/ds^2) kappa in the second data row

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;

% defining supported plate kernels

fkern1_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log',nu);           % build the desired kernel
fkern2_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth',nu);           % build the desired kernel

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';

% building system matrix

start = tic;
M = chunkermat(chnkr,fkern1_0, opts);
M2 = chunkermat(chnkr,fkern2_0, opts2);

c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + M2 + c0.*curv(:).^2.*eye(chnkr.npt);
M = M - 0.5*eye(2*chnkr.npt);

sys_0 = reshape(M,1,2*chnkr.npt,2*chnkr.npt);
sysmat1 = chunkermat(chnkr,fkern);
sysmat1 = reshape(sysmat1,nkappa,2*chnkr.npt,2*chnkr.npt);

sys = sysmat1 + sys_0;

t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)



%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa,d,sn,[],[],l,ising);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_bcs',kappa,d,sn,[],[],l,ising,nu);

skern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 's');

rhs = -bskern(src,chnkr);

% Solving linear system
sol = 0*rhs;
for i = 1:nkappa
sol(i:nkappa:end,:) = (squeeze(sys(i,:,:))\rhs(i:nkappa:end,:))*ws(i);
end

%%

uin = skern_0(src,targout_0);
if nkappa == 1
    uin = exp(1i*kappa(:).*nshift(iout).'*d).*skern(src,targout_0)*ws(1);
end

uscat = 0*uin;

ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_eval',kappa,d,sn,[],[],l,0,nu);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval',nu);

wts = repmat(chnkr.wts(:).',2,1);
nt = size(targout.r,2);

start1 = tic;
nbatch = ceil(2e5/chnkr.npt);
ntout = size(targout.r,2);

nshiftout = nshift(iout);
for i = 1:ceil(ntout/nbatch)
    iuse = ((i-1)*nbatch+1):min(ntout,i*nbatch);
    targi = []; targi.r = targout.r(:,iuse);
    nti = length(iuse);
    
    gevalmat_0 = chunkerkernevalmat(chnkr,ikern_0,targi,opts);
    gevalmat = ikern(chnkr,targi).* wts(:).';
    
    gevalmat = reshape(gevalmat,nkappa, nti, []);
    gevalmat = gevalmat + reshape(gevalmat_0,1,nti, []);
    gevalmat = exp(1i*kappa(:).*nshiftout(iuse)*d) .* gevalmat;
    % gevalmat = exp(1i*kappa(:).*nshift(iout)*d) .* gevalmat;
    
    gevalmat = reshape(permute(gevalmat, [2,1,3]), nti,[]);
    uscat(iuse,:) = gevalmat*sol;
end
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


utot = uscat+uin;



%%
chnkrs = [];
for i = -4:4
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = utot(:,1);

figure(1);clf
h = pcolor(X,Y, reshape(log10(abs(us)/norm(uin(:,1),inf)),size(X))); h.EdgeColor = 'None';
hold on
scatter(src.r(1,1),src.r(2,1),400,'r.')
plot(chnkrs,'k.','markersize',15)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');
% exportgraphics(gcf,'supported_acc.pdf','resolution',200)
% %%
figure(2);clf
us(iout) = utot(:,2);
% h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
h = pcolor(X,Y, reshape((imag(us)),size(X))); h.EdgeColor = 'None';
hold on
scatter(src.r(1,2),src.r(2,2),400,'r.')
plot(chnkrs,'k.','markersize',15)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');
% exportgraphics(gcf,'supported_sol.pdf','resolution',200)

