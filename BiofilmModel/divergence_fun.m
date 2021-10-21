function out_mat = divergence_fun(x_in_mat,y_in_mat)
x_part = par_x(x_in_mat);
y_part = par_y(y_in_mat);
out_mat = x_part + y_part;
end