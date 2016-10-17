load solution.dat;
load valX.dat;
load valY.dat;
N=size(valX,1);

[X,Y]=meshgrid(valX,valY)
mesh(X,Y,reshape(solution,N,N));
axis tight;



