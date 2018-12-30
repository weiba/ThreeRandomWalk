function auc=CaculateAUC(X,Y)
auc=0;
x(1)=1;y(1)=1;
row=size(X);
for i=2:row
auc=auc+(Y(i)+Y(i-1))*(X(i)-X(i-1))/2;
end;
auc=auc+Y(1)*X(1)/2;
end