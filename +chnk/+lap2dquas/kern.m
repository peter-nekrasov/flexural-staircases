function submat = kern(srcinfo,targinfo,type,kappa,d,s0,sn,l,ising,varargin)
%CHNK.LAP2D.KERN standard Laplace layer potential kernels in 2D
%
% Syntax: submat = kern(srcinfo,targinfo,type,varargin)

src = srcinfo.r(:,:);
targ = targinfo.r(:,:);

[~,ns] = size(src);
[~,nt] = size(targ);

switch lower(type)
% double layer
case {'d', 'double'}
    %srcnorm = chnk.normal2d(srcinfo);
    [~,grad] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    submat = -(grad(:,:,1).*srcinfo.n(1,:) + grad(:,:,2).*srcinfo.n(2,:));

% normal derivative of single layer
case {'sp', 'sprime'}
    targnorm = targinfo.n;
    [~,grad] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (grad(:,:,1).*nx + grad(:,:,2).*ny);

% Tangential derivative of single layer
case {'stau'}
    targnorm = targinfo.n;
    [~,grad] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    nx = repmat((targnorm(1,:)).',1,ns);
    ny = repmat((targnorm(2,:)).',1,ns);

    submat = (-grad(:,:,1).*ny + grad(:,:,2).*nx);

% Hilbert transform (two times the adjoint of stau)
case {'hilb'} 
    srcnorm = srcinfo.n;
    [~,grad] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    nx = repmat((srcnorm(1,:)),nt,1);
    ny = repmat((srcnorm(2,:)),nt,1);

    submat = 2*(grad(:,:,1).*ny - grad(:,:,2).*nx);

% tangential derivative Hilbert transform (tau of two times the adjoint of stau)
case {'hilbprime'} 
    srcnorm = srcinfo.n;
    targnorm = targinfo.n;
    [~,~,hess] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    nx = repmat((srcnorm(1,:)),nt,1);
    ny = repmat((srcnorm(2,:)),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);

    submat = 2*(-(hess(:,:,1).*ny - hess(:,:,2).*nx).*nytarg ...
        + (hess(:,:,2).*ny - hess(:,:,3).*nx).*nxtarg);

% normal derivative of double layer
case {'dprime','dp'}
    targnorm = targinfo.n;
    srcnorm = srcinfo.n;
    [~,~,hess] = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,ising);
    nxsrc = repmat(srcnorm(1,:),nt,1);
    nysrc = repmat(srcnorm(2,:),nt,1);
    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    submat = -(hess(:,:,1).*nxsrc.*nxtarg + hess(:,:,2).*(nysrc.*nxtarg+nxsrc.*nytarg)...
      + hess(:,:,3).*nysrc.*nytarg);

% single layer
case {'s', 'single'}
    submat = chnk.lap2dquas.green(src,targ,kappa,d,s0,sn,l,1);

otherwise
    error('Unknown quasiperiodic Laplace kernel type ''%s''.', type);
end

end

