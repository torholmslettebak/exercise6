nproc = [1,2,3,4,6];
runtime = [10.268501,5.148137,3.461355,2.609999,1.767766];

plot(nproc, runtime);
title ("MPI processes vs runtime");
xlabel ("nproc");
ylabel ("runtime[seconds]");
