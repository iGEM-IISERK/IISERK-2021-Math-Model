function out_mat = laplacian_fun(in_mat)
[x_comp,y_comp] = gradient_fun(in_mat);
out_mat = divergence_fun(x_comp,y_comp);
end
