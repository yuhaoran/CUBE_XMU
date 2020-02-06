function [ maxval,loc ] = maxloc( a )
[m,n]=size(a);
loc=[0,0];
maxval=-inf;
for j=1:n
  for i=1:m
    if a(i,j)>maxval
      maxval=a(i,j);
      loc=[i,j];
    end
  end
end

end

