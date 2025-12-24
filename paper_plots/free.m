%%
addpath(genpath('../../flexural-staircases'))
% zk = 0.8;
% zk = 0.3;
d = 1.2;
nu = 0.3; 

nnode = 62;
ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa = ts + amp*1i*sin(ts*d);
xip = 1 + amp*1i*d*cos(ts*d);
ws = ws*xip;

% kappa = kappa(1); ws = ws(1); ws = 1;
nkappa = length(kappa);


%%

nplot = 240;
nplot = 60;
xx = linspace(-6*d, 6*d,nplot);
yy = xx;
yy = linspace(0, 6*d,nplot/2) - 1.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
targmod = targ;
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);
%%

cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
cparams.maxchunklen = 2/zk;cparams.ifclosed = 1;cparams.eps = 1e-6;
nch = 20; A = 1;
% chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
% chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
chnkr = chunkerfunc(@(t) new_geom(t,d,A),cparams);
chnkr = reverse(chnkr);
% wtarg = cos_func(targmod.r(1,:),d,A) ;
wtarg = new_geom(targmod.r(1,:),d,A) ;
iout = targmod.r(2,:) > wtarg(2,:);

src = []; src.r = [[0;-2],[d/2;2]];
% src = []; src.r = [[0;-2],[d/2;1]];
% src.r = [0;-2];

targout = []; targout.r = targmod.r(:,iout);
targout_0 = []; targout_0.r = targ.r(:,iout);

%%

l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
[s0_l,sn_l] = chnk.lap2dquas.latticecoefs((1:N),d,kappa,l);
%%

ising = 0;
fkern1 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate',kappa,d,sn,s0_l,sn_l,l,ising,nu);
double = @(s,t) chnk.lap2dquas.kern(s,t,'d',kappa,d,s0_l,sn_l,l,ising);
hilbert = @(s,t) chnk.lap2dquas.kern(s,t,'hilb',kappa,d,s0_l,sn_l,l,ising);
opts = [];
opts.sing = 'smooth';

opts2 = [];
opts2.sing = 'smooth';

% building system matrix

start = tic;
sysmat1 = chunkermat(chnkr,fkern1, opts);
D = chunkermat(chnkr, double, opts);
H = chunkermat(chnkr, hilbert, opts2);     

sysmat1 = reshape(sysmat1,nkappa,4*chnkr.npt,2*chnkr.npt);
D = reshape(D,nkappa,chnkr.npt,chnkr.npt);
H = reshape(H,nkappa,chnkr.npt,chnkr.npt);


fkern1 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate',nu);
double = @(s,t) chnk.lap2d.kern(s,t,'d');
hilbert = @(s,t) chnk.lap2d.kern(s,t,'hilb');

opts = [];
opts.sing = 'log';

opts2 = [];
opts2.sing = 'pv';

% building system matrix

sysmat1_0 = chunkermat(chnkr,fkern1, opts);
D_0 = chunkermat(chnkr, double, opts);
H_0 = chunkermat(chnkr, hilbert, opts2); 

sysmat1_0 = reshape(sysmat1_0,1,4*chnkr.npt,2*chnkr.npt);
D_0 = reshape(D_0,1,chnkr.npt,chnkr.npt);
H_0 = reshape(H_0,1,chnkr.npt,chnkr.npt);

sysmat1 = sysmat1 + sysmat1_0; D = D + D_0; H = H + H_0;

D = permute(D,[2,3,1]);
H = permute(H,[2,3,1]);
s11b = permute(sysmat1(:,3:4:end,1:2:end),[2,3,1]);
s21b = permute(sysmat1(:,4:4:end,1:2:end),[2,3,1]);

k11tmp = permute(pagemtimes(s11b,H) -  2*((1+nu)/2)^2*pagemtimes(D,D),[3,1,2]);
k21tmp = permute(pagemtimes(s21b,H),[3,1,2]);

sysmat = zeros(nkappa,2*chnkr.npt,2*chnkr.npt);
sysmat(:,1:2:end,1:2:end) = sysmat1(:,1:4:end,1:2:end) + k11tmp;
sysmat(:,2:2:end,1:2:end) = sysmat1(:,2:4:end,1:2:end) + k21tmp;
% sysmat(:,1:2:end,1:2:end) = sysmat1(:,1:4:end,1:2:end) + sysmat1(:,3:4:end,1:2:end)*H  - 2*((1+nu)/2)^2*D*D;
% sysmat(:,2:2:end,1:2:end) = sysmat1(:,2:4:end,1:2:end) + sysmat1(:,4:4:end,1:2:end)*H;
sysmat(:,1:2:end,2:2:end) = sysmat1(:,1:4:end,2:2:end) + sysmat1(:,3:4:end,2:2:end);
sysmat(:,2:2:end,2:2:end) = sysmat1(:,2:4:end,2:2:end) + sysmat1(:,4:4:end,2:2:end);

D = [-1/2 + (1/8)*(1+nu).^2, 0; 0, 1/2];  % jump matrix 
D = reshape(kron(eye(chnkr.npt), D),1,2*chnkr.npt,2*chnkr.npt);

sys =  D + sysmat;
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

%%
ising = 1;
skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa,d,sn,s0_l,sn_l,l,ising);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa,d,sn,s0_l,sn_l,l,ising,nu);
skern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 's');

rhs = -bskern(src,chnkr);

% Solving linear system
sol = 0*rhs;
for i = 1:nkappa
sol(i:nkappa:end,:) = (squeeze(sys(i,:,:))\rhs(i:nkappa:end,:))*ws(i);
end
sol = reshape(sol,nkappa,[],size(sol,2));

dens_comb = zeros(nkappa,3*chnkr.npt,size(rhs,2));
dens_comb(:,1:3:end,:) = sol(:,1:2:end,:);
tmpsol = permute(sol,[2,3,1]);
dh = pagemtimes(H,tmpsol(1:2:end,:,:));
dh = permute(dh,[3,1,2]);
dens_comb(:,2:3:end,:) = dh;
dens_comb(:,3:3:end,:) = sol(:,2:2:end,:);
dens_comb = reshape(dens_comb,[],size(rhs,2));

%%

uin = skern_0(src,targout_0);
if nkappa == 1
    uin = exp(1i*kappa(:).*nshift(iout).'*d).*skern(src,targout_0)*ws(1);
end

uscat = 0*uin;

ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval',kappa,d,sn,s0_l,sn_l,l,0,nu);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval',nu);

wts = repmat(chnkr.wts(:).',3,1);
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
    uscat(iuse,:) = gevalmat*dens_comb;
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
% exportgraphics(gcf,'free_acc.pdf','resolution',200)
% %%
figure(2);clf
us(iout) = utot(:,2);
h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
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
% exportgraphics(gcf,'free_sol.pdf','resolution',200)


function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end