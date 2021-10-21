function out_mat = par_y(in_mat)
global xgridlen ygridlen;
out_mat = zeros(ygridlen+1,xgridlen+1);
for x = 0:xgridlen
    for y = 1:ygridlen-1
        lower_mid = (in_mat(y+1,x+1) + in_mat(y,x+1))/2;
        upper_mid = (in_mat(y+1,x+1) + in_mat(y+2,x+1))/2;
        out_mat(y+1,x+1) = (upper_mid-lower_mid);
    end
end
end