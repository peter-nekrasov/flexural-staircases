%%
addpath(genpath('../../flexural-staircases'))

% d = 1.2;
nu = 0.3; 
% zks = linspace(0.8,pi/d,40);
kappas = linspace(1,pi/d,10);
kappas = linspace(0.3+1e-2,pi/d,40);
% kappas = 1;
% poles = 0*kappas;
npoles = length(kappas);
poles = cell(1,npoles);



% cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
% cparams.maxchunklen = 2/(pi/d);cparams.ifclosed = 1;cparams.eps = 1e-6;
% nch = 20; A = 1;
% chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
% chnkr = reverse(chnkr);

tol = 1e-5;
for i = flip(1:npoles)
    i
kappa = kappas(i);
nkappa = 1;


kmin = 0.3; kmax = kappa-.1e-3;
poles{i} = neu_mode(chnkr,nu,kappa,kmin,kmax,d,tol);

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
for ii = 1:npoles
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

%%
kappa_vec = cell(3,1);
pole_vec = cell(3,1);
for i = 1:npoles
    for k = 1:size(pole_vec,1)
        if numel(poles{i}) > k-1
            kappa_vec{k} = [kappa_vec{k},kappas(i)];
            pole_vec{k} = [pole_vec{k},poles{i}(k)];
        end
    end
end


figure(2);clf

hold on
for ii = 1:size(pole_vec,1)
plot(kappa_vec{ii},real(pole_vec{ii}),'.-','linewidth',2,'MarkerSize',15)
end
plot(kappas,kappas,'k-','linewidth',1.5)
hold off
xlabel('$\xi_k$','interpreter','latex')
ylabel('$k$','interpreter','latex')
set(gca,'fontsize',18)
set(gca,'ticklabelinterpreter','latex')
drawnow()


function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end

function sys = neu_mat(chnkr,zk,nu,kappa,d)
% nkappa = length(kappa);
% 
% l=2; N = 40; a = 15; M = 1e4;
% sn = chnk.helm2dquas.latticecoefs((0:N).',zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
% %%
% quas_param = [];
% quas_param.kappa = kappa;
% quas_param.d = d;
% quas_param.l = l;
% quas_param.sn = sn;

start = tic;
% fkern1 =  @(s,t) chnk.helm2dquas.kern(zk, s, t, 'sprime',quas_param);
fkern1 = kernel('hq','sp',zk,kappa,d);
sysmat1 = chunkermat(chnkr,fkern1);

sys = -0.5*eye(size(sysmat1)) + sysmat1;
t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)


end

function modes = neu_mode(chnkr,nu,kappa,kmin,kmax,d,tol)


ncheb = 16;
xcheb = cos((2*(1:ncheb)-1)/2/ncheb*pi);
T = cos((0:(ncheb-1)).' .*acos(xcheb(:).')).';

% % get integrals of Chebyshev polynomials exactly
% [xleg,wleg] = lege.exps(ceil(ncheb/2)+3);
% Tleg = cos((0:(ncheb-1)).' .*acos(xleg(:).')).';
% wts = Tleg.'*wleg;


pan_ref_new = [kmin;kmax];
pans = [];
modes = [];

nref = 12;
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
            sys = neu_mat(chnkr,zks(i),nu,kappa,d);
            dets(i,1) = det(2*(sys));
            % sys = free_mat(chnkr,zks_l(i),nu,kappa,d);
            % dets(i,2) = det(2*squeeze(sys(1,:,:)));
            % sys = free_mat(chnkr,zks_r(i),nu,kappa,d);
            % dets(i,3) = det(2*squeeze(sys(1,:,:)));
        end

        coefs = (T\(dets.*sqrt(kappa-zks.').^0));
        % ints = wts(:).'*coefs;

        % check that the means of each half are comparable
        % idone =  max(abs(coefs(1,1))) <  min(max(abs(coefs(1,2:3))))*2;
        % idone = (abs(coefs(1,2)) <  2 *abs(coefs(1,3))) & (2*abs(coefs(1,2)) > abs(coefs(1,3))) ;

        % idone = idone & (max(abs(coefs(end-1:end,1)))/max(abs(coefs(1,1)))) < tol;
        idone = (max(abs(coefs(end-1:end,1)))/max(abs(coefs(1,1)))) < tol;
        % iref = abs(ints(1) - (ints(2)+ints(3))/2)/abs((ints(2)+ints(3))/2) < tol;
        % iref = abs(ints(1) - (ints(2)+ints(3))/2) < tol;
        if idone
            pans = [pans, pan_ref(:,k)];
            ipan_rm = [ipan_rm, k];

            c_cheb = T\(dets(:,1).*sqrt(kappa-zks.').^0);
            cs = c_cheb/c_cheb(end);
            
            B = .5*ones(ncheb-1,2);
            A = spdiags(B,[-1,1],ncheb-1,ncheb-1);
            A(1,2) = 1/sqrt(2);A(2,1) = 1/sqrt(2);
            en = zeros(1,ncheb-1); en(end)=1;
            cs(1) = sqrt(2)*cs(1);
            
            B = A - .5*cs(1:ncheb-1)*en;
            
            
            rts = eig(B);
            
            rts = rts(abs(rts)<1);
            rts = real(rts(abs(imag(rts))<1e-3));
            
            modes = [modes,diff(pan_ref(:,k))*(rts(:).'+1)/2+pan_ref(1,k)];
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

% toc;

% tic;
for k = 1:size(pan_ref,2)
    zks = diff(pan_ref(:,k))*(xcheb+1)/2+pan_ref(1,k);
    
    dets = zeros(ncheb,1);
    for i = 1:ncheb
        sys = neu_mat(chnkr,zks(i),nu,kappa,d);
        dets(i) = det(2*squeeze(sys(:,:)));
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
    rts = real(rts(abs(imag(rts))<1e-3));
    
    modes = [modes,diff(pan_ref(:,k))*(rts(:).'+1)/2+pan_ref(1,k)];
end
% toc;

modes = sort(modes);
end