close all
clear all
clc

del_t=100
end_time=1000

VidObj_live = VideoWriter("conc_profile_NisPV.avi");
open(VidObj_live);

for i = -10:0.1:0
    biocide_conc = 50*(10^i);
    fig = gen_conc_fig(biocide_conc, del_t, end_time)
    F1 = getframe(gcf);
    writeVideo(VidObj_live,F1);
    close(fig)
end
    
close(VidObj_live)
