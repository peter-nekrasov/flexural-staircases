%%
addpath(genpath('../../flexural-staircases'))

zks = linspace(0.8,pi/d,40);
zks = linspace(0.8,1.06,10);
poles = 0*zks;
npoles = length(zks);


cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
cparams.maxchunklen = 2/max(zks);cparams.ifclosed = 1;cparams.eps = 1e-6;
nch = 20; A = 1;
chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
chnkr = reverse(chnkr);


for i = 1:npoles

zk = zks(i);
d = 1.2;
nu = 0.3; 

amp = -0.3;

for j = 0:2
nleg = 32*2^j;

xs = cos((2*(1:nleg)-1)/2/nleg*pi);

tmin = zk+0.01; tmax = pi/d;
tmin = zk+.1; tmax = pi/d;
if i >1
tmin = max(tmin,poles(i-1));
end
tr = (tmax-tmin)*(xs+1)/2+tmin;
kappa = tr;

nkappa = length(kappa);

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

dets = zeros(nkappa,1);

for k = 1:nkappa
    dets(k) = det(2*squeeze(sys(k,:,:)));
end

%%
T = cos((0:(nleg-1)).' .*acos(xs(:).')).';

c_cheb = T\dets;

cs = c_cheb/c_cheb(end);

B = .5*ones(nleg-1,2);
A = spdiags(B,[-1,1],nleg-1,nleg-1);
A(1,2) = 1/sqrt(2);A(2,1) = 1/sqrt(2);
en = zeros(1,nleg-1); en(end)=1;
cs(1) = sqrt(2)*cs(1);

B = A - .5*cs(1:nleg-1)*en;


rts = eig(B);

rts = rts(abs(rts)<1);

rts= rts(abs(imag(rts))<1e-3);

rts = (tmax-tmin)*(rts+1)*.5 + tmin;

figure(5)
plot(kappa,abs(dets),'o-')

% figure(4)
% plot(rts,'o')
% title('Poles','Interpreter','latex')
% set(gca,'fontsize',16)

if abs(imag(rts)) < 1e-4
    poles(i) = rts;
    break
end
end
figure(1);clf
plot(zks(1:i),real(poles(1:i)),'.-','linewidth',2)
xlabel('$k$','interpreter','latex')
ylabel('$\xi_k$','interpreter','latex')
set(gca,'fontsize',18)
set(gca,'ticklabelinterpreter','latex')
drawnow()
end


%%
figure(1);clf
plot(zks,real(poles),'.-','linewidth',2)
xlabel('$k$','interpreter','latex')
ylabel('$\xi_k$','interpreter','latex')
set(gca,'fontsize',18)
set(gca,'ticklabelinterpreter','latex')

exportgraphics(gcf,'free_disp.pdf')




function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end