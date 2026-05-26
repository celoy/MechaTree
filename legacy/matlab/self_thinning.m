load SelfThinning.dat
 
gen = SelfThinning(:,1);
N   = SelfThinning(:,2);
M   = SelfThinning(:,3);

subplot(1,3,1)
hold on
plot(gen, N)
xlabel('generation')
ylabel('N')

subplot(1,3,2)
hold on
plot(gen, M)
xlabel('generation')
ylabel('M')

subplot(1,3,3)
hold on
loglog(N, M)
plot([10 1000],3*[100 .1],':r')
xlabel('N')
ylabel('M')
