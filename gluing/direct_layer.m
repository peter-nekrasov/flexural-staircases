function submat = direct_layer(src,targ,zk)
    [val,grad,hess,third] = chnk.flex2d.hkdiffgreen(zk,src.r(:,:),targ.r(:,:));  

    submat = zeros(size(targ.r(:,:),2),4*size(src.r(:,:),2));
    submat(:,1:4:end) = 1/(2*zk^2).*val;
    submat(:,2:4:end) = 1/(2*zk^2).*grad(:,:,1);
    submat(:,3:4:end) = 1/(2*zk^2).*hess(:,:,1);
    submat(:,4:4:end) = 1/(2*zk^2).*third(:,:,1);
end