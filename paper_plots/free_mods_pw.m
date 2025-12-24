%%
addpath(genpath('../../flexural-staircases'))

d = 1.2;
nu = 0.3; 
% zks = linspace(0.8,pi/d,40);
kappas = linspace(1,pi/d,10);
% kappas = 1;
% poles = 0*kappas;
npoles = length(kappas);
poles = cell(1,npoles);



cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
cparams.maxchunklen = 2/(pi/d);cparams.ifclosed = 1;cparams.eps = 1e-6;
nch = 20; A = 1;
chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
chnkr = reverse(chnkr);

tol = 1e-6;
for i = 10:npoles
kappa = kappas(i);
nkappa = 1;


kmin = 0.8; kmax = kappa-1e-2;
poles{i} = free_mode(chnkr,nu,kappa,kmin,kmax,d,tol);

figure(1);clf
hold on
for ii = 1:i
plot(kappas(ii)+0*real(poles{ii}),real(poles{ii}),'o','linewidth',2)
end
plot(kappas,kappas,'-','linewidth',2)
hold off
xlabel('$\xi_k$','interpreter','latex')
ylabel('$k$','interpreter','latex')
set(gca,'fontsize',18)
set(gca,'ticklabelinterpreter','latex')
drawnow()


end

%%
figure(2);clf
hold on
for ii = 1:i
plot(kappas(ii)+0*real(poles{ii}),real(poles{ii}),'o','linewidth',2)
end
plot(kappas,kappas,'-','linewidth',2)
hold off
xlabel('$\xi_k$','interpreter','latex')
ylabel('$k$','interpreter','latex')
set(gca,'fontsize',18)
set(gca,'ticklabelinterpreter','latex')
drawnow()
% exportgraphics(gcf,'free_disp.pdf')


function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end

function sys = free_mat(chnkr,zk,nu,kappa,d)
nkappa = length(kappa);

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


end

function modes = free_mode(chnkr,nu,kappa,kmin,kmax,d,tol)


ncheb = 16;
xcheb = cos((2*(1:ncheb)-1)/2/ncheb*pi);
T = cos((0:(ncheb-1)).' .*acos(xcheb(:).')).';

% % get integrals of Chebyshev polynomials exactly
% [xleg,wleg] = lege.exps(ceil(ncheb/2)+3);
% Tleg = cos((0:(ncheb-1)).' .*acos(xleg(:).')).';
% wts = Tleg.'*wleg;


pan_ref_new = [kmin;kmax];
pans = [];

nref = 10;
% tic;
for j = 1:nref
    pan_ref = pan_ref_new;
    ipan_rm = [];
    % tic;
    for k = 1:size(pan_ref,2)
        zks = diff(pan_ref(:,k))*(xcheb+1)/2+pan_ref(1,k);
        pan_ctr = mean(pan_ref(:,k));

        zks_l = diff(pan_ref(:,k))/2*(xcheb+1)/2+pan_ref(1,k);
        zks_r = diff(pan_ref(:,k))/2*(xcheb+1)/2+pan_ctr;

        dets = zeros(ncheb,3);

        for i = 1:ncheb
            sys = free_mat(chnkr,zks(i),nu,kappa,d);
            dets(i,1) = det(2*squeeze(sys(1,:,:)));
            sys = free_mat(chnkr,zks_l(i),nu,kappa,d);
            dets(i,2) = det(2*squeeze(sys(1,:,:)));
            sys = free_mat(chnkr,zks_r(i),nu,kappa,d);
            dets(i,3) = det(2*squeeze(sys(1,:,:)));
        end

        coefs = (T\(dets.*sqrt(kappa-zks.').^3));
        % ints = wts(:).'*coefs;

        % check that the means of each half are comparable
        idone =  max(abs(coefs(1,1))) <  min(max(abs(coefs(1,2:3))))*2;
        % idone = (abs(coefs(1,2)) <  2 *abs(coefs(1,3))) & (2*abs(coefs(1,2)) > abs(coefs(1,3))) ;

        idone = idone & (max(abs(coefs(end-1:end,1)))/max(abs(coefs(1,1)))) < tol;
        % iref = abs(ints(1) - (ints(2)+ints(3))/2)/abs((ints(2)+ints(3))/2) < tol;
        % iref = abs(ints(1) - (ints(2)+ints(3))/2) < tol;
        if idone
            pans = [pans, pan_ref(:,k)];
            ipan_rm = [ipan_rm, k];
        else
            a = pan_ref(1,k); b = pan_ref(2,k);
            pan_ref_new(:,k) = [a;pan_ctr];
            pan_ref_new = [pan_ref_new,[pan_ctr;b]];
        end

    end
    pan_ref_new = pan_ref_new(:, setdiff(1:size(pan_ref_new,2),ipan_rm));
    % toc;
end
pans = [pans,pan_ref];
modes = [];
% toc;

% tic;
for k = 1:size(pans,2)
    zks = diff(pans(:,k))*(xcheb+1)/2+pans(1,k);
    
    dets = zeros(ncheb,1);
    for i = 1:ncheb
        sys = free_mat(chnkr,zks(i),nu,kappa,d);
        dets(i) = det(2*squeeze(sys(1,:,:)));
    end

    c_cheb = T\(dets.*sqrt(kappa-zks.').^3);
    cs = c_cheb/c_cheb(end);
    
    B = .5*ones(ncheb-1,2);
    A = spdiags(B,[-1,1],ncheb-1,ncheb-1);
    A(1,2) = 1/sqrt(2);A(2,1) = 1/sqrt(2);
    en = zeros(1,ncheb-1); en(end)=1;
    cs(1) = sqrt(2)*cs(1);
    
    B = A - .5*cs(1:ncheb-1)*en;
    
    
    rts = eig(B);
    
    rts = rts(abs(rts)<1);
    rts = real(rts(abs(imag(rts))<1e-10));
    
    modes = [modes,diff(pans(:,k))*(rts+1)/2+pans(1,k)];
end
% toc;

modes = sort(modes);
end