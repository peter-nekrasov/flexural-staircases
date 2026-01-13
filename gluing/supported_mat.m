function [sys,sn,l] = supported_mat(chnkr,zk,nu,kappa,d)
nkappa = length(kappa);

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

end