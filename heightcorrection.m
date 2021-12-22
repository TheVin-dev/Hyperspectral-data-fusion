function corrected_heights = heightcorrection(lcmdata,show)
%% plane correctation  
heights = lcmdata.h_scaled;
[totalrows, totalcols] = size(heights);

lcmdata.showHeight()
[x,y] = ginput(3);
x = round(x,0);
y = round(y,0);

P = [x,y];
for point=1:length(P)
    h(point) = heights(P(point,2),P(point,1));
end

%%
tic
points = [P,h'];
syms x y z 
L = [x,y,z];
p1 = points(1,:);
p2 = points(2,:);
p3 = points(3,:);
normal = cross(p1-p2,p1-p3);
planeeq = dot(normal,L-p1);



zplane = solve(planeeq,z);
[x,y] = ndgrid(1:1:size(heights,2),1:1:size(heights,2));

zplaneFH = matlabFunction(zplane);
xlist = linspace(1,totalcols,totalcols);
ylist = xlist;
z = zplaneFH(x,y);

figure(1)
plot3(x,y,z2)
xlabel("x")
ylabel('y')
zlabel('height')
z = z(1:totalrows,1:totalcols);
% fm = fmesh(zplane,[1,totalcols, 1, totalrows]);
toc

if show

z = z * max(max(heights));
correct_heights = heights - z; 
figure(2)
subplot(1,3,1)
image(z,'CDatamapping','scaled');
daspect([1 1 1])
colorbar

subplot(1,3,2)
image(heights,'CDatamapping','scaled');
daspect([1 1 1])
colorbar

subplot(1,3,3)
image(correct_heights,'CDatamapping','scaled');
daspect([1 1 1])
colorbar

end
end