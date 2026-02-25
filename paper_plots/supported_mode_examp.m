zks = [1.180964293688935];
kappas = pi/d;
zks = 3.6;
kappas  =   2.682434068482174 - 2*pi/d;
nu = 0.3;

nplot = 200;
% nplot = 60;
xx = linspace(-4*d, 4*d,nplot);
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
figure(10);clf; t = tiledlayout('flow'); t.TileSpacing = 'tight';
figure(11);clf; t1 = tiledlayout('flow'); t1.TileSpacing = 'tight';


for j = 1:length(zks)

zk = real(zks(j)); kappa_rt = kappas(end); nkappa = 1;
if nkappa == 0, return, end
nkappa = length(kappa_rt);
chnkr = makedatarows(chnkr,2);
curv = signed_curvature(chnkr);
kp = arclengthder(chnkr,curv);
kpp = arclengthder(chnkr,kp);

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;


l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa_rt,(exp(1i*kappa_rt*d)),a,M,l+1);

ising = 0;
fkern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate',kappa_rt,d,sn,[],[],l,ising,nu);

opts2 = [];
opts2.quad = 'native';
opts2.sing = 'smooth';

% building system matrix

start = tic;
M3 = chunkermat(chnkr,fkern, opts2);   
M3 = reshape(M3,nkappa,2*chnkr.npt,2*chnkr.npt);

fkern_0l =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_log',nu);           % build the desired kernel
fkern_0s =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_smooth',nu);           % build the desired kernel

opts = [];
opts.sing = 'log';

M = chunkermat(chnkr,fkern_0l, opts);
M2 = chunkermat(chnkr,fkern_0s, opts2);

c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + M2 + c0.*curv(:).^2.*eye(chnkr.npt) - 0*eye(chnkr.npt); % extra term shows up for the general problem
M = M - 0.5*eye(2*chnkr.npt);
M = reshape(M,1,2*chnkr.npt,2*chnkr.npt);

sys = M3 + M;

t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

[u,sig,v] = svd(squeeze(sys));
dens = v(:,end);
sig(end,end)/sig(1,1)
%%

skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa_rt,d,sn,[],[],l,1);
% bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'free_plate_bcs',kappa_rt,d,sn,s0_l,sn_l,l,1,nu);

% Solving linear system
sol = dens;


ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_eval',kappa_rt,d,sn,[],[],l,0,nu);
ikern_0 = @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_eval',nu);

wts = repmat(chnkr.wts(:).',2,1);
nt = size(targout.r,2);

start1 = tic;
nbatch = ceil(2e5/chnkr.npt);
ntout = size(targout.r,2);
uscat = zeros(ntout,1);
nshiftout = nshift(iout);
for i = 1:ceil(ntout/nbatch)
    iuse = ((i-1)*nbatch+1):min(ntout,i*nbatch);
    targi = []; targi.r = targout.r(:,iuse);
    nti = length(iuse);
    
    gevalmat_0 = chunkerkernevalmat(chnkr,ikern_0,targi,opts);
    gevalmat = ikern(chnkr,targi).* wts(:).';
    
    gevalmat = reshape(gevalmat,nkappa, nti, []);
    gevalmat = gevalmat + reshape(gevalmat_0,1,nti, []);
    gevalmat = exp(1i*kappa_rt(:).*nshiftout(iuse)*d) .* gevalmat;
    % gevalmat = exp(1i*kappa(:).*nshift(iout)*d) .* gevalmat;
    
    gevalmat = reshape(permute(gevalmat, [2,1,3]), nti,[]);
    uscat(iuse,:) = gevalmat*sol;
end
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
plot(chnkrs,'k.','markersize',10)
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
plot(chnkrs,'k.','markersize',10)
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
