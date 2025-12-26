d = 1.2;
zk = 2;

cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
cparams.maxchunklen = 2/zk;cparams.ifclosed = 1;cparams.eps = 1e-6;
nch = 20; A = 1;
chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
chnkr = reverse(chnkr);

nplot = 40;
kr = linspace(1e-2,pi/d,nplot);
ki = linspace(-pi/d,pi/d,nplot);
[KR, KI] = meshgrid(kr,ki);
kappa = KR(:) + 1i*KI(:);

% kappa = linspace(1e-2,pi/d,20);

[sys, H, sn,s0_l, sn_l]= free_mat(chnkr,zk,nu,kappa,d);
targ = []; targ.r = [0;1.4];
src = targ;
uscat = free_rep(chnkr, zk, nu, kappa,d,sn,s0_l, sn_l,l, sys,H,src, targ);

%%
figure(1);clf
h = pcolor(KR,KI, reshape(imag(uscat), size(KR))); h.EdgeColor = 'none';

% figure(1);clf
% h = pcolor(KR,KI, reshape(angle(uscat), size(KR))); h.EdgeColor = 'none';



function [sys,H,sn,s0_l, sn_l,l] = free_mat(chnkr,zk,nu,kappa,d)
nkappa = length(kappa);
tic;
l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
[s0_l,sn_l] = chnk.lap2dquas.latticecoefs((1:N),d,kappa,l);
toc;
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


end


function uscat = free_rep(chnkr, zk,nu, kappa,d,sn,s0_l, sn_l,l, sys,H, src, targout)
% tic;
% l=2; N = 40; a = 15; M = 1e4;
% sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
% [s0_l,sn_l] = chnk.lap2dquas.latticecoefs((1:N),d,kappa,l);
% toc
nkappa = length(kappa);

ising = 1;
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa,d,sn,s0_l,sn_l,l,ising,nu);

rhs = -bskern(src,chnkr);

tic
sol = pagemldivide(permute(sys,[2,3,1]), permute(reshape(rhs,nkappa,[],size(rhs,2)), [2,3,1]));
sol = permute(sol,[3,1,2]); %sol = sol(:);
sol = reshape(sol,nkappa,[],size(rhs,2));
toc


dens_comb = zeros(nkappa,3*chnkr.npt,size(rhs,2));
dens_comb(:,1:3:end,:) = sol(:,1:2:end,:);
tmpsol = permute(sol,[2,3,1]);
dh = pagemtimes(H,tmpsol(1:2:end,:,:));
dh = permute(dh,[3,1,2]);
dens_comb(:,2:3:end,:) = dh;
dens_comb(:,3:3:end,:) = sol(:,2:2:end,:);
% dens_comb = reshape(dens_comb,[],size(rhs,2));
dens_comb = permute(dens_comb, [2,3,1]);

ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval',kappa,d,sn,s0_l,sn_l,l,0,nu);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval',nu);

wts = repmat(chnkr.wts(:).',3,1);
nt = size(targout.r,2);

start1 = tic;

    
gevalmat_0 = chunkerkernevalmat(chnkr,ikern_0,targout);
gevalmat = ikern(chnkr,targout).* wts(:).';

gevalmat = reshape(gevalmat,nkappa, nt, []);
gevalmat = gevalmat + reshape(gevalmat_0,1,nt, []);

gevalmat = permute(gevalmat, [2,3,1]);
uscat = pagemtimes(gevalmat,dens_comb);
uscat = uscat(:).';
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


end