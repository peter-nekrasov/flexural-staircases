function submat= kern(zk,srcinfo,targinfo,type,kappa,d,Sn,s0_l,sn_l,l,ising,varargin)

% 
% see also CHNK.FLEX2D.KERN
  
src = srcinfo.r;
targ = targinfo.r;

[~,ns] = size(src);
[~,nt] = size(targ);

%%% STANDARD LAYER POTENTIALS

switch lower(type)
case {'s', 'single'} % flexural wave single layer

   val = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
   submat = 1/(2*zk^2).*val;

case {'sp', 'sprime'} % normal derivative of flexural wave single layer

   targnorm = targinfo.n;
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);

   [~,grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l);  
   submat = 1/(2*zk^2).*(grad(:,:,1).*nxtarg + grad(:,:,2).*nytarg);

case {'d', 'double'} % normal derivative of flexural wave single layer

    srcnorm = srcinfo.n;
    nx = repmat((srcnorm(1,:)).',1,ns);
    ny = repmat((srcnorm(2,:)).',1,ns);
    
    [~,grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l);  
    submat = -1/(2*zk^2).*(grad(:,:,1).*nx + grad(:,:,2).*ny);


%%% CLAMPED PLATE KERNELS

% boundary conditions applied to a point source
case {'clamped_plate_bcs'}
    nxtarg = targinfo.n(1,:).'; 
    nytarg = targinfo.n(2,:).';  
    submat = zeros(2*nt,ns);
    
    [val, grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
    
    firstbc = 1/(2*zk^2).*val ;
    secondbc = 1/(2*zk^2).*(grad(:, :, 1).*nxtarg + grad(:, :, 2).*nytarg);
   
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

% kernels for the clamped plate integral equation
case {'clamped_plate'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);
   
   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);
   
   [~, ~, hess, third, fourth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
   
   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);
    
   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds;
   tauy = dy./ds;
   
   rx = targ(1,:).' - src(1,:);
   ry = targ(2,:).' - src(2,:);
   r2 = rx.^2 + ry.^2;

   rn = rx.*nx + ry.*ny;
   rtau = rx.*taux + ry.*tauy;
   ntargtau = nxtarg.*taux + nytarg.*tauy;

   rntarg = rx.*nxtarg + ry.*nytarg;

   K11 = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny))) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

   K12 = -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

   K21 = -(1/(2*zk^2).*(fourth(:, :, 1).*(nx.*nx.*nx.*nxtarg) + fourth(:, :, 2).*(nx.*nx.*nx.*nytarg + 3*nx.*nx.*ny.*nxtarg) + ...
          fourth(:, :, 3).*(3*nx.*nx.*ny.*nytarg + 3*nx.*ny.*ny.*nxtarg) + fourth(:, :, 4).*(3*nx.*ny.*ny.*nytarg +ny.*ny.*ny.*nxtarg)+...
          fourth(:, :, 5).*(ny.*ny.*ny.*nytarg)) ) - ...
          (3/(2*zk^2).*(fourth(:, :, 1).*(nx.*taux.*taux.*nxtarg)+ fourth(:, :, 2).*(nx.*taux.*taux.*nytarg + 2*nx.*taux.*tauy.*nxtarg + ny.*taux.*taux.*nxtarg) +...
          fourth(:, :, 3).*(2*nx.*taux.*tauy.*nytarg + ny.*taux.*taux.*nytarg + nx.*tauy.*tauy.*nxtarg + 2*ny.*taux.*tauy.*nxtarg) + ...
          fourth(:, :, 4).*(nx.*tauy.*tauy.*nytarg +2*ny.*taux.*tauy.*nytarg + ny.*tauy.*tauy.*nxtarg) +...
          fourth(:, :, 5).*(ny.*tauy.*tauy.*nytarg)));

   K22 = -(1/(2*zk^2).*(third(:,:, 1).*(nx.*nx.*nxtarg) +third(:, :, 2).*(nx.*nx.*nytarg + 2*nx.*ny.*nxtarg) + third(:, :, 3).*(2*nx.*ny.*nytarg + ny.*ny.*nxtarg)+...
         third(:, :,4).*(ny.*ny.*nytarg))) + ...
         (1/(2*zk^2).*(third(:,:, 1).*(taux.*taux.*nxtarg) +third(:, :, 2).*(taux.*taux.*nytarg + 2*taux.*tauy.*nxtarg) + third(:, :, 3).*(2*taux.*tauy.*nytarg + tauy.*tauy.*nxtarg)+...
         third(:, :,4).*(tauy.*tauy.*nytarg)));

  submat = zeros(2*nt,2*ns);
  
  submat(1:2:end,1:2:end) = K11;
  submat(1:2:end,2:2:end) = K12;
    
  submat(2:2:end,1:2:end) = K21;
  submat(2:2:end,2:2:end) = K22;

% clamped plate kernels for plotting
case {'clamped_plate_eval'}

    submat = zeros(nt,2*ns);

    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [~, ~, hess, third] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  

    K1 = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;

% clamped plate kernels for far field evaluation
case {'clamped_plate_eval_ff'}

    submat = zeros(nt,2*ns);

    srcnorm = srcinfo.n;
    srctang = srcinfo.d;
    nx = repmat(srcnorm(1,:),nt,1);
    ny = repmat(srcnorm(2,:),nt,1);
    dx = repmat(srctang(1,:),nt,1);
    dy = repmat(srctang(2,:),nt,1);
    ds = sqrt(dx.*dx+dy.*dy);

    taux = dx./ds;
    tauy = dy./ds;

    [~, ~, hess, third] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  

    K1 = -(1/(2*zk^2).*(third(:, :, 1).*(nx.*nx.*nx) + third(:, :, 2).*(3*nx.*nx.*ny) +...
       third(:, :, 3).*(3*nx.*ny.*ny) + third(:, :, 4).*(ny.*ny.*ny)) ) - ...
       (3/(2*zk^2).*(third(:, :, 1).*(nx.*taux.*taux) + third(:, :, 2).*(2*nx.*taux.*tauy + ny.*taux.*taux) +...
       third(:, :, 3).*(nx.*tauy.*tauy + 2*ny.*taux.*tauy) + third(:, :, 4).*(ny.*tauy.*tauy)));  % G_{ny ny ny} + 3G_{ny tauy tauy}

    K2 =  -(1/(2*zk^2).*(hess(:, :, 1).*(nx.*nx) + hess(:, :, 2).*(2*nx.*ny) + hess(:, :, 3).*(ny.*ny)))+...
          (1/(2*zk^2).*(hess(:, :, 1).*(taux.*taux) + hess(:, :, 2).*(2*taux.*tauy) + hess(:, :, 3).*(tauy.*tauy))); % -G_{ny ny}  + G_{tauy tauy}

    submat(:,1:2:end) = K1;
    submat(:,2:2:end) = K2;

%%% FREE PLATE KERNELS

% boundary conditions applied to a point source
case {'free_plate_bcs'}
    targnorm = targinfo.n;
    targtang = targinfo.d;
    targd2 = targinfo.d2;
    nu = varargin{1};
    
    [~, ~, hess, third] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  

    nxtarg = repmat((targnorm(1,:)).',1,ns);
    nytarg = repmat((targnorm(2,:)).',1,ns);
    
    dx1 = repmat((targtang(1,:)).',1,ns);
    dy1 = repmat((targtang(2,:)).',1,ns);
    
    ds1 = sqrt(dx1.*dx1+dy1.*dy1); 
    
    d2x1 = repmat((targd2(1,:)).',1,ns);
    d2y1 = repmat((targd2(2,:)).',1,ns);
    
    tauxtarg = dx1./ds1;
    tauytarg = dy1./ds1;
    
    denom = sqrt(dx1.^2 + dy1.^2).^3;
    numer = dx1.*d2y1 - d2x1.*dy1;
    
    kappatarg = numer ./ denom; % target curvature
    
    firstbc = 1/(2*zk^2).*(hess(:, :, 1).*(nxtarg.*nxtarg) + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*(nytarg.*nytarg))+...
    nu/(2*zk^2).*(hess(:, :, 1).*(tauxtarg.*tauxtarg) + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*(tauytarg.*tauytarg));
    
    secondbc = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
    third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))+...
    (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
    third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg+ tauytarg.*tauytarg.*nxtarg) +...
    + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg))+...
    (1-nu).*kappatarg.*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
    (1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)));

    submat = zeros(2*nt,ns);
    submat(1:2:end,:) = firstbc;
    submat(2:2:end,:) = secondbc;

% kernels for the free plate integral equation 
case {'free_plate'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   targnorm = targinfo.n;
   targtang = targinfo.d;
   targd2 = targinfo.d2;
   nu = varargin{1};

   [~, ~, hess, third, fourth] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  

   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   nxtarg = repmat((targnorm(1,:)).',1,ns);
   nytarg = repmat((targnorm(2,:)).',1,ns);


   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy); 

   taux = dx ./ ds;
   tauy = dy ./ ds;

   dx1 = repmat((targtang(1,:)).',1,ns);
   dy1 = repmat((targtang(2,:)).',1,ns);

   ds1 = sqrt(dx1.*dx1+dy1.*dy1); 

   d2x1 = repmat((targd2(1,:)).',1,ns);
   d2y1 = repmat((targd2(2,:)).',1,ns);

   tauxtarg = dx1./ds1;
   tauytarg = dy1./ds1;

   denom = sqrt(dx1.^2 + dy1.^2).^3;
   numer = dx1.*d2y1 - d2x1.*dy1;

   kappatarg = numer ./ denom; % target curvature

   hilb = chnk.lap2dquas.kern(srcinfo,targinfo,'hilb',kappa,d,s0_l,sn_l,l,0);
   hilbp = chnk.lap2dquas.kern(srcinfo,targinfo,'hilbprime',kappa,d,s0_l,sn_l,l,0);
   hilb = (1+nu) * hilb / 2;
   hilbp = (1+nu) * hilbp / 2;
   
   K11 = -(1/(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
        third(:, :, 4).*(nytarg.*nytarg.*ny))) - ...
       nu./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ;  % first kernel with no hilbert transforms (G_{nx nx ny + nu G_{taux taux ny}).

   K12 =  1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg)+...
           nu/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg) ;    % G_{nx nx} + nu G_{taux taux}
   
   K21 = kappatarg./(2*zk^2).*(1-nu).*((third(:, :, 1).*(nxtarg.*nxtarg.*nx) + third(:, :, 2).*(nxtarg.*nxtarg.*ny + 2*nxtarg.*nytarg.*nx) +...
            third(:, :, 3).*(2*nxtarg.*nytarg.*ny + nytarg.*nytarg.*nx) +...
            third(:, :, 4).*(nytarg.*nytarg.*ny)) - ...
           (third(:, :, 1).*(tauxtarg.*tauxtarg.*nx) + third(:, :, 2).*(tauxtarg.*tauxtarg.*ny + 2*tauxtarg.*tauytarg.*nx) +...
            third(:, :, 3).*(2*tauxtarg.*tauytarg.*ny + tauytarg.*tauytarg.*nx) +...
            third(:, :, 4).*(tauytarg.*tauytarg.*ny)) ) ...
        - 1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*ny + 3*nxtarg.*nxtarg.*nytarg.*nx) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*ny + 3*nxtarg.*nytarg.*nytarg.*nx) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*ny +nytarg.*nytarg.*nytarg.*nx)+...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*ny)) - ...
          ((2-nu)/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*nx) + fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*ny + tauxtarg.*tauxtarg.*nytarg.*nx + 2*tauxtarg.*tauytarg.*nxtarg.*nx)+...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*ny + 2*tauxtarg.*tauytarg.*nxtarg.*ny + tauytarg.*tauytarg.*nxtarg.*nx + 2*tauxtarg.*tauytarg.*nytarg.*nx)+...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*ny + 2*tauxtarg.*tauytarg.*nytarg.*ny + tauytarg.*tauytarg.*nytarg.*nx) +...
          fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*ny)) ) ...          
          -hilbp/2;

   K22 = 1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*nxtarg) + third(:, :, 2).*(3*nxtarg.*nxtarg.*nytarg) +...
       third(:, :, 3).*(3*nxtarg.*nytarg.*nytarg) + third(:, :, 4).*(nytarg.*nytarg.*nytarg))  +...
        (2-nu)/(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg) + third(:, :, 2).*(tauxtarg.*tauxtarg.*nytarg + 2*tauxtarg.*tauytarg.*nxtarg) +...
        third(:, :, 3).*(2*tauxtarg.*tauytarg.*nytarg + tauytarg.*tauytarg.*nxtarg) +...
        + third(:, :, 4).*(tauytarg.*tauytarg.*nytarg)) + ... % G_{nx nx nx} + (2-nu) G_{taux taux nx}
        + kappatarg.*(1-nu).*(1/(2*zk^2).*(hess(:, :, 1).*tauxtarg.*tauxtarg + hess(:, :, 2).*(2*tauxtarg.*tauytarg) + hess(:, :, 3).*tauytarg.*tauytarg)-...
           1/(2*zk^2).*(hess(:, :, 1).*nxtarg.*nxtarg + hess(:, :, 2).*(2*nxtarg.*nytarg) + hess(:, :, 3).*nytarg.*nytarg));
    
   K11H =  -(1+ nu)/2*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)) )  - ...
       (1+ nu)/2*nu.*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy))) + (1+ nu)/4*hilb;

   K21H = kappatarg.*(1-nu).*(-((1+ nu)/2).*(1./(2*zk^2).*(third(:, :, 1).*(tauxtarg.*tauxtarg.*taux) + third(:, :, 2).*(tauxtarg.*tauxtarg.*tauy + 2*tauxtarg.*tauytarg.*taux) +...
       third(:, :, 3).*(2*tauxtarg.*tauytarg.*tauy + tauytarg.*tauytarg.*taux) +...
        third(:, :, 4).*(tauytarg.*tauytarg.*tauy)))  + ...
        ((1+ nu)/2)*(1./(2*zk^2).*(third(:, :, 1).*(nxtarg.*nxtarg.*taux) + third(:, :, 2).*(nxtarg.*nxtarg.*tauy+ 2*nxtarg.*nytarg.*taux) +...
        third(:, :, 3).*(2*nxtarg.*nytarg.*tauy +nytarg.*nytarg.*taux) +...
        third(:, :, 4).*(nytarg.*nytarg.*tauy)))) ...
        -(1+ nu)/2.*(1/(2*zk^2).*(fourth(:, :, 1).*(nxtarg.*nxtarg.*nxtarg.*taux) + ...
          fourth(:, :, 2).*(nxtarg.*nxtarg.*nxtarg.*tauy + 3*nxtarg.*nxtarg.*nytarg.*taux) + ...
          fourth(:, :, 3).*(3*nxtarg.*nxtarg.*nytarg.*tauy + 3*nxtarg.*nytarg.*nytarg.*taux) +...
          fourth(:, :, 4).*(3*nxtarg.*nytarg.*nytarg.*tauy + nytarg.*nytarg.*nytarg.*taux) +...
          fourth(:, :, 5).*(nytarg.*nytarg.*nytarg.*tauy)) ) - ...
          ((2-nu)/2)*(1+nu).*(1/(2*zk^2).*(fourth(:, :, 1).*(tauxtarg.*tauxtarg.*nxtarg.*taux) + ...
          fourth(:, :, 2).*(tauxtarg.*tauxtarg.*nxtarg.*tauy + tauxtarg.*tauxtarg.*nytarg.*taux + 2*tauxtarg.*tauytarg.*nxtarg.*taux) + ...
          fourth(:, :, 3).*(tauxtarg.*tauxtarg.*nytarg.*tauy + 2*tauxtarg.*tauytarg.*nxtarg.*tauy + tauytarg.*tauytarg.*nxtarg.*taux + 2*tauxtarg.*tauytarg.*nytarg.*taux) +...
          fourth(:, :, 4).*(tauytarg.*tauytarg.*nxtarg.*tauy + 2*tauxtarg.*tauytarg.*nytarg.*tauy + tauytarg.*tauytarg.*nytarg.*taux) +...
         fourth(:, :, 5).*(tauytarg.*tauytarg.*nytarg.*tauy))) ;
    
    submat = zeros(4*nt,2*ns);
    
    submat(1:4:end,1:2:end) = K11;
    submat(1:4:end,2:2:end) = K12;
    
    submat(2:4:end,1:2:end) = K21;
    submat(2:4:end,2:2:end) = K22;
    
    submat(3:4:end,1:2:end) = K11H;
    submat(4:4:end,1:2:end) = K21H;

% free plate kernels used for plotting 
case {'free_plate_eval'}
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   nu = varargin{1};

   [val,grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds; 
   tauy = dy./ds;
  
   K1 = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 
   K1H = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*(taux) + grad(:, :, 2).*tauy));                    % G_{tauy}
   K2 = 1/(2*zk^2).*val;

   submat = zeros(nt,3*ns);
   submat(:,1:3:end) = K1;
   submat(:,2:3:end) = K1H;
   submat(:,3:3:end) = K2;

% free plate kernels for far field evaluation
case {'free_plate_eval_ff'}           
   srcnorm = srcinfo.n;
   srctang = srcinfo.d;
   nu = varargin{1};

   [val,grad] = chnk.flex2dquas.green(src,targ,zk,kappa,d,Sn,l,0);  
   nx = repmat(srcnorm(1,:),nt,1);
   ny = repmat(srcnorm(2,:),nt,1);

   dx = repmat(srctang(1,:),nt,1);
   dy = repmat(srctang(2,:),nt,1);

   ds = sqrt(dx.*dx+dy.*dy);

   taux = dx./ds; 
   tauy = dy./ds;
  
   K1 = (-1/(2*zk^2).*(grad(:, :, 1).*(nx) + grad(:, :, 2).*ny)); 
   K1H = ((1 + nu)/2).*(-1/(2*zk^2).*(grad(:, :, 1).*(taux) + grad(:, :, 2).*tauy));                    % G_{tauy}
   K2 = 1/(2*zk^2).*val;

   submat = zeros(nt,3*ns);
   submat(:,1:3:end) = K1;
   submat(:,2:3:end) = K1H;
   submat(:,3:3:end) = K2;

end

if ising == 1
    submat = submat + chnk.flex2d.kern(zk,srcinfo,targinfo,type,varargin{:});
end
end


