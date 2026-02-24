d = 2;
zk = 7;

cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
cparams.eps = 1e-10;
cparams.maxchunklen = 2/zk;cparams.ifclosed = 1;
nch = 20; A = -0.5;
% chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
chnkr0 = chunkerfunc(@(t) cos_func(t,d,A),cparams);
chnkr1 = reverse(chnkr0);

cparams = rmfield(cparams,'ta');
cparams = rmfield(cparams,'tb');
chnkr2 = chunkerfunc(@(t) starfish(t,3),cparams); 
chnkr2 = move(chnkr2,[], [0;2.5],0.3,0.5);
src = []; src.r = [0;0]; src.n = [1;0];



rend = chunkends(chnkr,[1,chnkr1.nch]);
rend = rend(:,[2,3]);

verts = [rend, rend - [0;10]];
edge2verts = [[1;2], [4;3], [3;1], [2;4],[NaN;NaN]];
fchnk = cell(1,5);
fchnk{1} = chnkr0;
fchnk{5} = chnkr2;
cgrph = chunkgraph(verts,edge2verts,fchnk);

cparams.ifclosed = 0;
 cparams.ta = -d/2; cparams.tb = d/2;
chnkr0 = chunkerfunc(@(t) cos_func(t,d,A),cparams);

chnkr = merge([chnkr1, chnkr2]);
chnkrplot = merge([chnkr0, chnkr2]);
figure(1);clf
plot(chnkrplot,'linewidth',2)
axis equal
chnkr.npt

nu = 0.3; 

nnode = 62;
ts = linspace(-pi/d,pi/d,nnode);
ts = ts(2:end);
ws = 1/(nnode-1);

amp = -0.3;
kappa = ts + amp*1i*sin(ts*d);
xip = 1 + amp*1i*d*cos(ts*d);
ws = ws*xip;

% kappa = pi/d + 1e-1; ws = 1;
% kappa = kappa(10); ws = ws(1); ws = 1;
nkappa = length(kappa);

chnkr = makedatarows(chnkr,2);
curv = signed_curvature(chnkr);
kp = arclengthder(chnkr,curv);
kpp = arclengthder(chnkr,kp);

chnkr.data(1,:,:) = kp;
chnkr.data(2,:,:) = kpp;

%%

nplot = 240;
% nplot = 60;
xx = linspace(-4*d, 4*d,nplot);
yy = xx;
yy = linspace(0, 6*d,3*nplot/3) - 1.2;
% yy = linspace(0, 6*d,3*nplot/4) - 3.2;
[X,Y] = meshgrid(xx,yy);
targ = []; targ.r = [X(:).'; Y(:).'];

targmod = [];
targmod.r = real([mod(targ.r(1,:)+d/2,d)-d/2;targ.r(2,:)]);
% targmod = targ;
nshift = round((targ.r(1,:)-targmod.r(1,:))/d);
%%

% cparams = []; cparams.ta = -d/2; cparams.tb = d/2;
% cparams.maxchunklen = 2/zk;cparams.ifclosed = 1;cparams.eps = 1e-6;
% nch = 20; A = 1;
% % chnkr = chunkerfuncuni(@(t) cos_func(t,d,A),nch,cparams);
% chnkr = chunkerfunc(@(t) cos_func(t,d,A),cparams);
% chnkr = chunkerfunc(@(t) new_geom(t,d,A),cparams);
% cparams.ta = -1; cparams.tb = 1; cparams.ifclosed = 0;
% rs = [(-1:0.3:-0.1).' 0*(-1:0.3:-0.1).'+1; -0.1 0.8; -0.5 0.6; -2/3 -1/2; 0 -1; 2/3 -1/2; 2/3 1/2; 1/2 0.8; (0.5:0.1:1).' (0.5:0.1:1).'*0+1;].';
% rs(1,:) = d/2*rs(1,:);
% coefs = get_splines(rs);
% chnkr = chunkerfunc(@(t) geom_eval(t,coefs),cparams);
% chnkr = reverse(chnkr);
% wtarg = cos_func(targmod.r(1,:),d,A) ;
% wtarg = new_geom(targmod.r(1,:),d,A) ;
% iout = targmod.r(2,:) > wtarg(2,:);

src = []; src.r = [[0;-2],[d/2;1.5]];
% src = []; src.r = [[0;-2],[d/2;1]];
% src.r = [0;-2];

chnkrs = [];
for i = (-6:6)
    chnkrs = [chnkrs, chnkrplot + [i*d;0]];
end
chnkrs = merge(chnkrs);
% iout = ~chunkerinterior(chnkrs,targ);
iout = chunkgraphinregion(cgrph,targmod)==1;

targout = []; targout.r = targmod.r(:,iout);
targout_0 = []; targout_0.r = targ.r(:,iout);

% figure(1); clf 
% plot(chnkr)
% hold on
% quiver(chnkr)
% drawnow

%%

l=2; N = 40; a = 15; M = 1e4;
sn = chnk.flex2dquas.latticecoefs((0:N).',zk,d,kappa,(exp(1i*kappa*d)),a,M,l+1);
[s0_l,sn_l] = chnk.lap2dquas.latticecoefs((1:N),d,kappa,l);
%

alpha = reshape(exp(1i*kappa(:)*d),[nkappa,1,1]); 
nsub = 1;

ising = 0;
fkern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate',kappa,d,sn,[],[],l,ising,nu,nsub);

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

c0 = (nu - 1)*(nu + 3)*(2*nu - 1)/(2*(3 - nu));

M = chunkermat(chnkr,fkern_0l, opts);
M2 = chunkermat(chnkr,fkern_0s, opts2);

M = M - 0.5*eye(2*chnkr.npt);
M(2:2:end,1:2:end) = M(2:2:end,1:2:end) + c0.*curv(:).^2.*eye(chnkr.npt) - 0*eye(chnkr.npt); % extra term shows up for the general problem

M = reshape(M,1,2*chnkr.npt,2*chnkr.npt);
M2 = reshape(M2,1,chnkr.npt,chnkr.npt);

for ii = -nsub:nsub
    if ii ~= 0 
        Msub = chunkerkernevalmat(chnkr + ii*[d;0],fkern_0l,chnkr);
        M2sub = chunkerkernevalmat(chnkr + ii*[d;0],fkern_0s,chnkr,struct('forcesmooth',true));
  
        Msub = reshape(Msub,[1,2*chnkr.npt,2*chnkr.npt]);
        M2sub = reshape(M2sub,[1,chnkr.npt,chnkr.npt]);

        M = M + alpha.^(ii).*Msub;
        M2 = M2 + alpha.^(ii).*M2sub;
    end
end

M(:,2:2:end,1:2:end) = M(:,2:2:end,1:2:end) + M2 ; % extra term shows up for the general problem

sys = M3 + M;

t1 = toc(start);
fprintf('%5.2e s : time to assemble matrix\n',t1)

%%
ising = 1;
skern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 's',kappa,d,sn,s0_l,sn_l,l,ising);
bskern =  @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_bcs',kappa,d,sn,s0_l,sn_l,l,ising,nu);
skern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 's');
bskern_0 =  @(s,t) chnk.flex2d.kern(zk, s, t, 'supported_plate_bcs',nu);

rhs = -bskern(src,chnkr);

% Solving linear system
sol = 0*rhs;
for i = 1:nkappa
sol(i:nkappa:end,:) = (squeeze(sys(i,:,:))\rhs(i:nkappa:end,:))*ws(i);
end

%%

uin = skern_0(src,targout_0);
if nkappa == 1
    uin = exp(1i*kappa(:).*nshift(iout).'*d).*skern(src,targout)*ws(1);
end

uscat = 0*uin;

ikern = @(s,t) chnk.flex2dquas.kern(zk, s, t, 'supported_plate_eval',kappa,d,sn,s0_l,sn_l,l,0,nu);
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

us = (NaN+NaN*1i)*zeros(1,size(targ.r,2));
us(iout) = utot(:,1);

f1=figure(1);clf
f1.Position = [1 1 643 441];
C = reshape(log10(abs(us)/norm(uin(:,1),inf)), size(X)); 
h = pcolor(X,Y,C);
h.EdgeColor = 'None'; 
h.FaceColor = 'texturemap'; 
h.AlphaData = ~isnan(C);
h.FaceAlpha = 'texturemap';
hold on
scatter(src.r(1,1),src.r(2,1),400,'r.')
plot(chnkrs,'k-','LineWidth',2.5)
c = colorbar;
% clim([1e-6 1e-5])
% c.Ticks = floor(c.Limits(1)) : ceil(c.Limits(2));
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
vv = sort(abs(C(:)));
clim([min(C(:)),-vv(3)])
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');
exportgraphics(f1,'supported_acc.pdf','ContentType','vector','Resolution',300);
% %%
f2=figure(2);clf
f2.Position = [1 1 643 441];
us(iout) = utot(:,2);
C = reshape((imag(us)),size(X));
h = pcolor(X,Y,C);
h.EdgeColor = 'None'; 
h.FaceColor = 'texturemap'; 
h.AlphaData = ~isnan(C);
h.FaceAlpha = 'texturemap';
hold on
scatter(src.r(1,2),src.r(2,2),300,'r.')
plot(chnkrs,'k-','LineWidth',2.5)
c = colorbar;
hold off
axis equal
xlim([min(X(:)),max(X(:))])
ylim([min(Y(:)),max(Y(:))])
set(gca,'FontSize',16)
set(gca,'TickLabelInterpreter','latex');
set(c,'TickLabelInterpreter','latex');
exportgraphics(f2,'supported_sol.pdf','ContentType','vector','Resolution',300);

load('gong.mat')
sound(y)

function [r,d,d2] = cos_func(t,d,A)
% parameterization of sinusoidal boundary with period d and amplitude A
omega = 2*pi/d;
r = [t(:), A*cos(omega*t(:))].';
d = [ones(length(t),1), -omega*A*sin(omega*t(:))].';
d2 = [zeros(length(t),1), -omega^2*A*cos(omega*t(:))].';
end
