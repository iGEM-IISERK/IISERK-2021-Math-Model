antimicrobe_con_d = 50e-3
del_t= 100
end_time = [200, 400, 600, 800];

for t = end_time
    filename = string(t)+"dead.png"
    fig = gen_conc_fig(antimicrobe_con_d,del_t,t)
    saveas(gcf, filename)
    close
    
end

