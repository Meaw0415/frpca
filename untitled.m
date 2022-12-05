errf_svd_L = [6.284328665055552e-08,6.284328665055552e-08,6.284328665055552e-08,6.284328665055552e-08,6.284328665055552e-08];
errf_svd_S = [1.704397447849009e-10,1.704397447849009e-10,1.704397447849009e-10,1.704397447849009e-10,1.704397447849009e-10];
figure;
subplot(1, 3, 1)
plot(1:my_ticks, errf_svd_L, '*-', 'LineWidth', 2);
hold on;
plot(1:my_ticks, errf_svd_S, '+-', 'LineWidth', 2);
hold off;
title("L AND S ERROR RPCA VS FRPCA");
legend("fast", "svd");
xticks(1:my_ticks);
xticklabels(num2cell(T));
xlabel("T");
ylabel("error (mse)");