function out_mat = par_x(in_mat)
global xgridlen ygridlen;
out_mat = zeros(ygridlen+1,xgridlen+1);
for y = 0:ygridlen
    for x = 1:xgridlen-1
        lower_mid = (in_mat(y+1,x+1) + in_mat(y+1,x))/2;
        upper_mid = (in_mat(y+1,x+1) + in_mat(y+1,x+2))/2;
        out_mat(y+1,x+1) = (upper_mid-lower_mid);
    end
end
end