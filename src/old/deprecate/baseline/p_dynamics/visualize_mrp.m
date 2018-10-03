close all
figure(); 
hold on

test_vector = [1;0;0];

N = 100;
[X,Y,Z] = sphere(N);
CO(:,:,1) = zeros(N); % red
CO(:,:,2) = zeros(N);  % green
CO(:,:,3) = zeros(N); % blue
surf(X,Y,Z,CO)

traces = zeros(3, size(MP_state,1));
for i = 1:size(MP_state,1)
    R = mrp2dcm(MP_state(i,1:3));
    traces(:,i) = R*test_vector;
end
plot3(traces(1,:), traces(2,:), traces(3,:), 'w');

xi = mrp2dcm(test_state(1:3))*test_vector;  scatter3(xi(1), xi(2), xi(3), 'gd');
xf = mrp2dcm(x_eq(1:3))*test_vector;        scatter3(xf(1), xf(2), xf(3), 'rx');
   
%grid on
%axis equal
%set(findall(gcf,'type','text'),'FontSize',32);set(gca,'FontSize',32)
