close all
% plot the projection of 4 dimensional parameter to ks and kd
clear;
% y = load('/Users/cmhshirley/Desktop/hill_qnstop/HQ_likely_region_Aconcentration/sOutput.txt');
% y = load('/Users/cmhshirley/Desktop/hill_qnstop/HQ_norm_region/sOutput.txt');
y = load('./20190529_fig/sim.txt');

alpha = 0.3;
beta = 0.7;
gamma = 0.3;

num = 20;
niter = 100;

size(y)

p = 4;
ite = y(:,2);
f_xi = y(:, 3);
f_bt = y(:, 4);
f_mx = y(:, 5);
xi = y(:, 6:5+p);
bt = y(:, 6+p:5+2*p);
TAU = y(:, 6+2*p);
W = y(:,7+2*p:6+2*p+p*p);
stability = y(:, end);
% Q = y(:,7+2*p+p*p:6+2*p+2*p*p);
% Eig = y(:, 7+2*p+2*p*p:end);

clf;
figure(1)
ite = 1:1:niter*num;
plot(ite, f_xi, 'r-', ite, f_bt, 'b-', ite, f_mx, 'k-')
xlabel('Iteration')
ylabel('f(X)', 'Interpreter','latex')
legend('Ellipsoid Center','Best Sampled', 'Worst Sampled')
% axis([0 100 0.73 1])

figure(2)
% single
% oneit = 1:1:niter;
% plot(oneit, f_xi(oneit), 'r-',  ...
%      oneit, f_bt(oneit), 'b--',  ...
%      oneit, f_mx(oneit), 'k:',  'LineWidth', 2)

%  avg
avg_f_xi = reshape(f_xi,[niter,num]);
avg_f_bt = reshape(f_bt,[niter,num]);
avg_f_mx = reshape(f_mx,[niter,num]);
plot(1:niter, mean(avg_f_xi, 2), 'r-', ...
     1:niter,mean(avg_f_bt, 2), 'b--', ...
     1:niter,mean(avg_f_mx, 2), 'k:', 'LineWidth', 2)
% plot(oneit, f_xi(oneit), 'r-', oneit, f_bt(oneit), 'b-', oneit, f_mx(oneit), 'k-')
xlabel('Iteration')
ylabel('$f(X)$', 'Interpreter','latex')
legend('Ellipsoid Center','Best Sampled', 'Worst Sampled')

% figure(3)
% plot(ite,stability, 'k-')
% xlabel('Iteration')
% ylabel('Stability', 'Interpreter','latex')

figure(4)
avg_sta = reshape(stability,[niter,num]);
plot(1:niter,mean(avg_sta, 2), 'k-')
xlabel('Iteration')
ylabel('Stability', 'Interpreter','latex')
legend('Max Probability-sim')

fmin = min(f_bt);

figure(5)
% for j = 0:num-1
for j = 0:num-1
%     fmin = min(f_bt(j*niter+1:j*niter+niter));
%     fmax = max(f_bt(j*niter+1:j*niter+niter));
    for i = 1: niter
       if stability(i+j*niter)>=beta && f_bt(i+j*niter)<=(1+gamma)*fmin 
            hold on
            pause(0.01)
            
            ww = reshape(W(i+j*niter,:), [p,p]);
            J = ww(1:2,1:2);
            LL= ww(1:2,3:4);
            L = ww(3:4,1:2);
            K = ww(3:4,3:4);
            Kinv = inv(K);

            ftest = @(a,b) [a-xi(i+j*niter,1), b-xi(i+j*niter,2)]*...
                (J - LL*Kinv*L)*...
                [a-xi(i+j*niter,1); b-xi(i+j*niter,2)]-TAU(i+j*niter)^2;
            ez3 = ezplot(ftest, [-3 3 -3 3]);
            set(ez3,'color',[0 0 1])
%             scatter(bt(i+j*niter,1),bt(i+j*niter,2))
%             scatter(xi(i+j*niter,1),xi(i+j*niter,2))
            i
        end      
    end
j
end

box on;
title('')
xlabel('$\log_{10}(k_a)$', 'Interpreter','latex')
ylabel('$\log_{10}(k_d)$', 'Interpreter','latex')
axis([-3 3 -3 3])

