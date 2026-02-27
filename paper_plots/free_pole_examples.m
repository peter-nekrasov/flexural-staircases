% load('free_disp.mat')
% load('free_disp2.mat')
% load('square_poles1.mat')
nplot = 200;
% nplot = 60;
xx = linspace(-6*d, 6*d,nplot);
yy = xx;
yy = linspace(0, 4*d,nplot/2) - 1.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
% targmod = targ;
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);

% wtarg = cos_func(targmod.r(1,:),d,A) ;
% iout = targmod.r(2,:) > wtarg(2,:);
iout = chunkgraphinregion(cgrph,targmod)==1;
targout = []; targout.r = targmod.r(:,iout);
targout_0 = []; targout_0.r = targ.r(:,iout);


%%
pole_vec = [];
zks_vec = [];
for i = 1:length(poles)
    if ~isempty(poles{i})
        pole_vec = [pole_vec, poles{i}(1)];
        % zks_vec = [zks_vec, zks(i) + 0* poles{i}];
    else
        pole_vec = [pole_vec, NaN];
        % zks_vec = [zks_vec, zks(i)];
    end
    zks_vec = [zks_vec, zks(i)];
end

inan = find(isnan(pole_vec));

ibdry = zeros(2,length(inan)+1);
ibdry(:,1) = [1; inan(1) - 1];
for i = 2:length(inan)
    ibdry(:,i) = [inan(i-1) + 1; inan(i)-1];
end
ibdry(:,end) = [inan(end)+1; length(pole_vec)];

ibdry(:,ibdry(1,:)>ibdry(2,:)) = [];

% iexamp = round(mean(ibdry,1));
iexamp = ibdry(2,:) - 2;
% % iexamp = [iexamp([1,3]),89,90];
% iexamp = [iexamp([1,3,4,5])];

% iexamp = 1:length(pole_vec);
% iexamp = iexamp(~isnan(pole_vec));

figure(5);clf
plot(zks_vec, real(pole_vec),'o-','linewidth',2)
hold on
% for ii = 1:length(zks)
% plot(zks(ii)+0*real(poles{ii}),real(poles{ii}),'o','linewidth',2)
% end
plot(zks,zks,'-','linewidth',2)
hold off
xlabel('$k$','interpreter','latex')
ylabel('$\xi_k$','interpreter','latex')
set(gca,'fontsize',18)
set(gca,'ticklabelinterpreter','latex')
drawnow()
% exportgraphics(gcf,'free_disp.pdf')

%%
figure(10);clf; t = tiledlayout('flow'); t.TileSpacing = 'tight';
figure(11);clf; t1 = tiledlayout('flow'); t1.TileSpacing = 'tight';


for j = 1:length(iexamp)

zk = zks_vec(iexamp(j)); kappa_rt = pole_vec(iexamp(j)); nkappa = 1;
if nkappa == 0, return, end
l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa_rt,(exp(1i*kappa_rt*d)),a,M,l+1);
[s0_l,sn_l] = chnk.lap2dquas.latticecoefs((1:N),d,kappa_rt,l);


ising = 0;
fkern1 =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate',kappa_rt,d,sn,s0_l,sn_l,l,ising,nu);
double = @(s,t) chnk.lap2dquas.kern(s,t,'d',kappa_rt,d,s0_l,sn_l,l,ising);
hilbert = @(s,t) chnk.lap2dquas.kern(s,t,'hilb',kappa_rt,d,s0_l,sn_l,l,ising);
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

[u,sig,v] = svd(squeeze(sys));
dens = v(:,end);
sig(end,end)/sig(1,1)
%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa_rt,d,sn,s0_l,sn_l,l,1);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa_rt,d,sn,s0_l,sn_l,l,1,nu);

% Solving linear system
sol = dens;


ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_eval',kappa_rt,d,sn,s0_l,sn_l,l,0,nu);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'free_plate_eval',nu);


dens_comb = zeros(3*chnkr.npt,1);
dens_comb(1:3:end) = sol(1:2:end);
dens_comb(2:3:end) = H*sol(1:2:end);
dens_comb(3:3:end) = sol(2:2:end);

wts = repmat(chnkr.wts(:).',3,1);

start1 = tic;
uscat = chunkerkerneval(chnkr, ikern_0,dens_comb,targout);
uscat = uscat + ikern(chnkr,targout) * (dens_comb .* wts(:));
uscat = uscat.*exp(1i*kappa_rt*d*nshift(iout).');
t2 = toc(start1);
fprintf('%5.2e s : time for kernel eval (for plotting)\n',t2)


%%
chnkrs = [];
for i = -10:10
    chnkrs = [chnkrs, chnkr + [i*d;0]];
end
chnkrs = merge(chnkrs);

figure(10);nexttile()
us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = uscat;
h = pcolor(X,Y, reshape((real(us)),size(X))); h.EdgeColor = 'None';
hold on
plot(chnkrs,'k.','markersize',15)
c = colorbar;
c.Label.String = '$\Im v_\xi$';
c.Label.Interpreter = 'latex';
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');

figure(11);nexttile()
h = pcolor(X,Y, reshape((abs(us)),size(X))); h.EdgeColor = 'None';
hold on
plot(chnkrs,'k.','markersize',15)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
% xlabel('$x_1$','Interpreter','latex')
% ylabel('$x_2$','Interpreter','latex')
c.Label.String = '$|v_\xi|$';
c.Label.Interpreter = 'latex';
set(gca,'FontSize',18)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');

end
