
rs = [(-1:0.3:-0.1).' 0*(-1:0.3:-0.1).'+1; -0.1 0.8; -0.5 0.6; -2/3 -1/2; 0 -1; 2/3 -1/2; 2/3 1/2; 1/2 0.8; (0.5:0.1:1).' (0.5:0.1:1).'*0+1;].';
% rs = [3.1113    3.4245    3.6456    3.2587    3.4429    4.5483    4.3457    4.6220    5.0089;...
%     1.0105    0.9552    0.7342    0.1630   -0.8318   -0.1686    0.6789    0.9368    1.0289];
% figure(5)
% a = drawpolygon;
% rs = a.Position.';
% % rs = rs(:,1:end-1);
% 
% rs = [5.1536    5.4113    5.5402    5.5770    5.5586    5.5034    5.6322    6.3870    6.6263     6.4607    6.3502    6.4054    6.5711    6.8657;...
%     1.5078    1.5262    1.4894    1.2869    0.8819    0.2375   -0.1491   -0.0938    0.1271    0.4768    0.9371    1.3421    1.4894    1.5078];
rs = [-1.1221     -0.62345      -.9482     -0.73426     -0.90047    -0.20946     -0.34643     -0.01401      0.24454      0.38147      0.74317      0.39228;...
       1.4804       1.2218      0.35383     -0.10786     -0.9117      -0.8727      -0.1448      0.24303     -0.27407     -0.82811     -0.16327      0.64932];

% chnkr = sort(chnkr);
% rs = chnkr.r(:,:);
rs = rs(:,1:end-1);
ts = linspace(-1,1,length(rs)+1); ts = ts(1:end-1);

ncheb = floor((size(rs,2)-1)/2);
% ncheb =7;


% tsub = ts;
A4 = [1+0*ts.',cos(pi*ts.'*(1:ncheb)), sin(pi*ts.'*(1:ncheb))];

cs4 = A4\(rs - ts.*[1;0]).';
rcheb = (A4*cs4).' + ts.*[1;0];


figure(3);clf
plot(rs(1,:), rs(2,:),'x')
hold on
plot(rcheb(1,:), rcheb(2,:),'.')
hold off
axis equal


rend = FP2([-1,1],cs4);
d = diff(rend(1,:));

figure(4);clf
plot(rcheb(1,:), rcheb(2,:),'.')
hold on
plot(rcheb(1,:)+d, rcheb(2,:),'.')
plot(rcheb(1,:)+2*d, rcheb(2,:),'.')
hold off


cparams = []; cparams.ta = -1; cparams.tb = 1; cparams.maxchunklen = 1/(pi/d);
chnkr = chunkerfunc(@(t) FP2(t,cs4),cparams);


cparams.ifclosed = 0;
chnkr0 = chunkerfunc(@(t) FP2(t,cs4),cparams);

figure(5);clf
plot(chnkr,'.')
hold on
for i = 1:3
plot(chnkr + i*[d;0],'.')
end
hold off
axis equal
chnkr.npt
chnkr = chnkr + [-mean(rend(1,:));0];
rend = rend - [mean(rend(1,:));0];

verts = [rend, rend - [0;10]];
edge2verts = [[1;2], [4;3], [3;1], [2;4]];
fchnk = cell(1,4);
fchnk{1} = chnkr0;
cgrph = chunkgraph(verts,edge2verts,fchnk);


chnkr = reverse(chnkr);
figure(6);clf
plot(cgrph,'.')
hold on
plot(chnkr,'x')
hold off
axis equal

% function [r,d,d2] = FP2(ts, cs)
function [r] = FP2(ts, cs)
ncheb = (size(cs,1))/2;
r  = ([1+0*ts(:),cos(pi*ts(:)*(1:ncheb)), sin(pi*ts(:)*(1:ncheb))]*cs).' + ts(:).'.*[1;0];
d  = ([0*ts(:),-sin(pi*ts(:)*(1:ncheb)).*(1:ncheb)*pi, cos(pi*ts(:)*(1:ncheb)).*(1:ncheb)*pi]*cs).'+ [1;0];
d2 = ([0*ts(:),-cos(pi*ts(:)*(1:ncheb)).*(1:ncheb).^2*pi^2, -sin(pi*ts(:)*(1:ncheb)).*(1:ncheb).^2*pi^2]*cs).';

end

