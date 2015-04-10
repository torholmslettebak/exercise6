problemsize = 2048;
threads = [1,2,3,4,6];

runtime = [10.244505, 5.691899, 4.184082, 3.388541, 2.621776];

plot(threads, runtime);
title ("number of threads vs runtime for n = 2048");
xlabel ("threads");
ylabel ("runtime[seconds]");
