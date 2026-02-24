
figure(5);clf
a = drawpolygon;
rs = a.Position.';
if mod(length(rs),2) == 0; rs = rs(:,1:end-1); end

%% person
% rs = 1.4*[0.33813      0.10541      0.18145      0.38422      0.37039      0.20449      0.4      0.55012      0.73906      0.84044...
%     0.66993      0.64228       0.8174        0.875      0.65      0.53477      0.58007      0.47638      0.39343;...
%       0.80365      0.59051      0.47664      0.63139      0.30438     0.091241     0.07       0.2     0.047445      0.12628...
%       0.31606      0.59927      0.44453      0.56131      0.80365      0.83365      0.95      1.1       0.9492];
% rs = rs-mean(rs,2)+[1;0];
% rs =[ 1.5333      0.81271       0.5785       1.8486       3.5782       3.8754       3.9024       3.1548       2.4071;...
%       0.85112      0.33767       1.4096       2.4906       2.3825       1.6168      0.17552      0.49981       1.5898];

% 
% %% fish
% rs = [1.4733       1.2888       1.0783      0.4      0.4      1.0696 1.2714       1.4646       1.45       1.45;...
%       0.87834      0.68301      0.68084      0.8      0.27716      0.4    0.4      0.31623      0.5     0.7];
% rs(2,:) = rs(2,:) - mean(rs(2,:));

%% blob
% rs = [1.5       1.0678      0.61818      0.58758       1.0584       1.5       1.4;...
%       0.86916      0.66201      0.92566      0.27596      0.53725      0.34894      0.61964];
% rs = rs(:,1:end-1);
% 

% chnkr = sort(chnkr);
% rs = chnkr.r(:,:);
% rs = rs(:,1:end-1);
ts = linspace(-1,1,length(rs)+1); ts = ts(1:end-1);

ncheb = floor((size(rs,2)-1)/2);
% ncheb =7;


% tsub = ts;
A4 = [1+0*ts.',cos(pi*ts.'*(1:ncheb)), sin(pi*ts.'*(1:ncheb))];

cs4 = A4\(rs).';
rcheb = (A4*cs4).';

%%

rcheb = FP2(linspace(-1,1,100), cs4);

figure(3);clf
plot(rs(1,:), rs(2,:),'x')
hold on
plot(rcheb(1,:), rcheb(2,:),'-')
hold off
axis equal

%%
rend = FP2([-1,1],cs4);
% d = diff(rend(1,:));
d = 2;

figure(4);clf
plot(rcheb(1,:), rcheb(2,:),'.')
hold on
plot(rcheb(1,:)+d, rcheb(2,:),'.')
plot(rcheb(1,:)-d, rcheb(2,:),'.')
hold off
axis equal

cparams = []; cparams.ta = -1; cparams.tb = 1; cparams.maxchunklen = 1/(pi/d);
cparams.eps = 1e-8;
chnkr = chunkerfunc(@(t) FP2(t,cs4),cparams);


% figure(5);clf
% plot(chnkr,'.')
% hold on
% for i = 1:3
% plot(chnkr + i*[d;0],'.')
% end
% hold off
% axis equal
chnkr.npt

% cgrph = tochunkgraph(chnkr);

figure(6);clf
% plot(cgrph,'.')
% hold on
plot(chnkr,'x')
% hold off
axis equal

% figure(7);clf
% plot_regions(cgrph)

% function [r,d,d2] = FP2(ts, cs)
function [r] = FP2(ts, cs)
ncheb = (size(cs,1))/2;
r  = ([1+0*ts(:),cos(pi*ts(:)*(1:ncheb)), sin(pi*ts(:)*(1:ncheb))]*cs).';
d  = ([0*ts(:),-sin(pi*ts(:)*(1:ncheb)).*(1:ncheb)*pi, cos(pi*ts(:)*(1:ncheb)).*(1:ncheb)*pi]*cs).';
d2 = ([0*ts(:),-cos(pi*ts(:)*(1:ncheb)).*(1:ncheb).^2*pi^2, -sin(pi*ts(:)*(1:ncheb)).*(1:ncheb).^2*pi^2]*cs).';

end

