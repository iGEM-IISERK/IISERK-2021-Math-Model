function conc_fig = gen_conc_fig(antimicrobe_con_d,del_t,end_time)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

max = 1; min = 0;
bacPerCell = 100;
time_arr=[];
mean_live =[]; mean_bact = []; mean_dead =[];

% Length Dimensionality: 1mm = 128 grid pixels (aka xgridlen)
global xgridlen ygridlen del_t;

%% Characteristic length
h0 = 1*10^(-3);

%% Grid dimensions
xgridlen = 128;
ygridlen = 128;

%% differential lengths
del_x = h0/xgridlen;
del_y = h0/ygridlen;

%% number of iteration
loop_count = end_time/del_t;
% Seting up Constants
%% Outer hemispherical dome radius (i.e. dome of 0.1 bacterical fraction)
set_radius = 38;
biofilm_layers = 4;
large_num = 10000;

%% Densities of Solveny and Bacteria
den_solv = 1000;
den_bact = 1000;

%% Gravity
gravity = 10;

%%Aqueous diffusion coefficient of Nisin, D0
nisin_diffusion_coeff_D0 = 0.05*10^(-9);%NisinA
% nisin_diffusion_coeff_D0 = 0.32*10^(-9);%NisinPV


%% Maximum biocide consumption rate for nisin
% k1 = 2*10^(-7); %NisinA
k1 = 14*10^(-7); %NisinPV

%%Biocide Monod half saturation coefficient
d0 = 1*10^(-3);

%% Mobility parameter (capital lambda)
lamb = 1*10^(-9);

%% Maximum decay rate of live bacteria due to biocide action
% k2 = 2*10^(-5); %NisinA
k2 = 14*10^(-5); %Nisin PV


%% Live bacteria decay Monod half-saturation coefficient
kd = 1*10^(-5);

%% Boltzmann constant
kb = 1.38065*10^(-23);

%% Temperature
temp_T = 296;

%% Distortional energy
tau1 = 8*10^(6);

%% Strength of the bulk mixing energy
tau2 = 3*10^(17);

%% Generalized polymerization parameter
cap_N = 1*10^(3);

%% Flory-Huggins mixing parameter
chi = 0.55;

%% Viscosity of the bacteria and the solvent
nu_bac = 5;
nu_sol = 1.002*10^(-3);

%% % Setting up the initial Grid Map
%% Initial grid set up
gridmap = csvread("0grid.csv");
bact_vol_frac = gridmap/bacPerCell;
for y = 1:set_radius+1
    for x = 0:xgridlen
        if(gridmap(y+1,x+1) == 1 || gridmap(y+1,x+1) == 2 || gridmap(y+1,x+1) == 3 || gridmap(y+1,x+1) == 4)
            bact_vol_frac(y+1,x+1) = -0.1*sqrt((x-(xgridlen/2))*(x-(xgridlen/2)) + y*y)+3.9;
        end
    end
end

sol_vol_frac = 1 - bact_vol_frac;
bact_live_vol_frac = bact_vol_frac;
bact_dead_vol_frac = bact_vol_frac - bact_live_vol_frac;

in_b = bact_vol_frac;
in_l = bact_live_vol_frac;
in_d = bact_dead_vol_frac;
% Setting up the initial velocity matrices
% Maximum speed at the left and right top points in the flow is 2cm/s ( = 2x10x128 pixels/s = 2560pixels/s)
%%Maxminum fluid speed
max_speed = 2*10^(-2);

%%Velocity Matries
x_vel = zeros(ygridlen+1,xgridlen+1);
y_vel = zeros(ygridlen+1,xgridlen+1);

%%Boundary conditions on velocity:
%% left and right
for y = 0:ygridlen
    x_vel(y+1,1) = (max_speed/(ygridlen*ygridlen))*y*y;
    x_vel(y+1,xgridlen+1) = (max_speed/(ygridlen*ygridlen))*y*y;
end

%% USING CONSERVATION OF VOLUME FLOW RATE
cen_max_speed = max_speed*ygridlen/(ygridlen-set_radius+3);

%% central
for y = set_radius-4:ygridlen
    x_vel(y+1,(xgridlen/2)+1) = (cen_max_speed/((ygridlen-set_radius+4)*(ygridlen-set_radius+4)))*(y-set_radius+4)*(y-set_radius+4);
end

%% Bi-quadratic X-component velocity variation in NON Film Region
for y = set_radius-3:ygridlen
    he = x_vel(y+1,1);
    hc = x_vel(y+1,(xgridlen/2)+1);
    for x = 1:(xgridlen/2)-1
        x_vel(y+1,x+1) = ((hc-he)*(((16/xgridlen^4)*((x-(xgridlen/2))^4)) - ((8/xgridlen^2)*((x-(xgridlen/2))^2)) +1)) + he;
    end
end

%% Bi-quadratic X-component velocity variation in Film Region
for y = 1:set_radius-4
    he = x_vel(y+1,1);
    for x = 1:xgridlen
        if(gridmap(y+1,x+1) == 10)
            len = 2*x;
            break
        end
    end
    for x =1:len/2
        x_vel(y+1,x+1) = ((0-he)*(((16/len^4)*((x-(len/2))^4)) - ((8/len^2)*((x-(len/2))^2)) +1)) + he;
    end
end

%% Y-component velocity of layer of 0.4 bacterial fraction (i.e. inner most layer)
for y = 1:set_radius
    for x = 1:(xgridlen/2)-1
        if(gridmap(y+1,x+1) == 4)
            slope = -1*(x-xgridlen/2)/y;
            y_vel(y+1,x+1) = slope*x_vel(y+1,x+1);
        end
    end
end

%% Y-component velocity obeying divergence condition for rest of the biofilm layers
for y = 1:set_radius+1
    for x = 1:(xgridlen/2)-1
        if(gridmap(y+1,x+1) == 1 || gridmap(y+1,x+1) == 2 || gridmap(y+1,x+1) == 3)
            y_vel(y+1,x+1) = y_vel(y,x+1) - x_vel(y+1,x+1) + x_vel(y+1,x);
        end
    end
end

%% Quadratic Y-component velocity variation in Film Region
for x = (xgridlen/2) - set_radius:(xgridlen/2)-1
    first = 1;
    for y = 1:ygridlen
        if(gridmap(y+1,x+1) == 0 && first == 1)
            first = 0;
            maxvy = y_vel(y,x+1);
            ht = y-1;
        end
        if(first == 0)
            y_vel(y+1,x+1) = (-1*maxvy/((ht - ygridlen)^2))*(y-ygridlen)*(y-2*ht+ygridlen);
        end
    end
end

%% Quadratic Y-component velocity variation in NON Film Region
ht = (xgridlen/2) - set_radius;
for y = 1:ygridlen-1
    maxvy = y_vel(y+1,ht+1);
    for x = 1:ht
        y_vel(y+1,x+1) = (-1*maxvy/ht^2)*x*(x-2*ht);
    end
end

%% Setting up velocities in right half
for y = 1:ygridlen
    for x = 1:(xgridlen/2)-1
        x_vel(y+1,129-x) = x_vel(y+1,1+x);
        y_vel(y+1,129-x) = -1*y_vel(y+1,1+x);
    end
end
%quiver(x_vel,y_vel)
% Main program loop


video_speed_up_rate = 5;
video_length = end_time/video_speed_up_rate;
total_frames = 24*video_length;
frame_every_x_loops = round(loop_count/total_frames);


for counter = 0:loop_count
    %% Setting up all the matrices
    
    %% Density Matrix
    density_mat = den_bact*bact_vol_frac + den_solv*sol_vol_frac;
    
    %% Pressure Matrix
    pressure_mat = zeros(ygridlen+1,xgridlen+1);
    for x = 0:xgridlen
        for y = 0:ygridlen
            speed_sq = x_vel(y+1,x+1)*x_vel(y+1,x+1) + y_vel(y+1,x+1)*y_vel(y+1,x+1);
            pressure_mat(y+1,x+1) = -1*density_mat(y+1,x+1)*(speed_sq + gravity*y/(1000*ygridlen));
        end
    end
    
    %% Effective Diffusion Coefficient (De) Matrix
    diffusion_coeff_mat_De = zeros(ygridlen+1,xgridlen+1);
    for x = 0:xgridlen
        for y = 0:ygridlen
            diffusion_coeff_mat_De(y+1,x+1) = (2*nisin_diffusion_coeff_D0*(1-0.4*bact_vol_frac(y+1,x+1)))/((2+0.4*bact_vol_frac(y+1,x+1))*(1+29*bact_vol_frac(y+1,x+1)));
        end
    end
    
    %% Antimicrobial Concentration Matrix
    antimicrobe_con_mat = zeros(ygridlen+1,xgridlen+1);
    for x = 0:xgridlen
        for y = 0:ygridlen
            if(sol_vol_frac(y+1,x+1) ~= 0)
                antimicrobe_con_mat(y+1,x+1) = antimicrobe_con_d;
            end
        end
    end
    
    %% Viscosity Matrix
    viscosity_mat = nu_sol*sol_vol_frac + nu_bac*bact_vol_frac;
    
    %% Chemical Potential df/d(phi_b) Matrix
    chem_pot_mat = zeros(ygridlen+1,xgridlen+1);
    lap_bact_vol_frac = laplacian_fun(bact_vol_frac);
    
    for x = 0:xgridlen
        for y = 0:ygridlen
            if(bact_vol_frac(y+1,x+1) ~= 0 && bact_vol_frac(y+1,x+1) ~= 1)
                chem_pot_mat(y+1,x+1) = kb*temp_T*((-1*tau1*lap_bact_vol_frac(y+1,x+1))+(tau2*(((1+log(bact_vol_frac(y+1,x+1)))/cap_N)-log(1-bact_vol_frac(y+1,x+1))-1+(chi*(1-2*bact_vol_frac(y+1,x+1))))));
            end
        end
    end
    
    % All Partials of needed variables
    vel_x_par_x = par_x(x_vel);
    vel_x_par_y = par_y(x_vel);
    vel_y_par_x = par_x(y_vel);
    vel_y_par_y = par_y(y_vel);
    [press_grad_x, press_grad_y] = gradient_fun(pressure_mat);
    [bact_frac_grad_x, bact_frac_grad_y] = gradient_fun(bact_vol_frac);
    [con_grad_x, con_grad_y] = gradient_fun(antimicrobe_con_mat);
    [chem_grad_x, chem_grad_y] = gradient_fun(chem_pot_mat);
    
    %% Solvent fraction change equation
    sol_con = sol_vol_frac.*antimicrobe_con_mat;
    
    x_term_1 = sol_vol_frac.*(x_vel.*antimicrobe_con_mat - diffusion_coeff_mat_De.*con_grad_x);
    y_term_1 = sol_vol_frac.*(y_vel.*antimicrobe_con_mat - diffusion_coeff_mat_De.*con_grad_y);
    
    con_par_t = -1*divergence_fun(x_term_1,y_term_1) - k1*(bact_live_vol_frac.*antimicrobe_con_mat)./(d0+antimicrobe_con_mat);
    
    sol_con = sol_con + del_t*con_par_t;
    
    %% Live bacteria fraction change equation
    live_par_t = -1*divergence_fun(bact_live_vol_frac.*x_vel, bact_live_vol_frac.*y_vel) + lamb*(divergence_fun(bact_live_vol_frac.*chem_grad_x, bact_live_vol_frac.*chem_grad_y)) - k2*(bact_live_vol_frac.*antimicrobe_con_mat)./(kd+antimicrobe_con_mat);
    
    bact_live_vol_frac = bact_live_vol_frac + del_t*live_par_t;
    
    %% Dead bacteria fraction change equation
    dead_par_t = -1*divergence_fun(bact_dead_vol_frac.*x_vel, bact_dead_vol_frac.*y_vel) + lamb*(divergence_fun(bact_dead_vol_frac.*chem_grad_x, bact_dead_vol_frac.*chem_grad_y)) + k2*(bact_live_vol_frac.*antimicrobe_con_mat)./(kd+antimicrobe_con_mat);
    
    bact_dead_vol_frac = bact_dead_vol_frac + del_t*dead_par_t;
    
    %% Strain Tensor Matrix matrices
    strain_tensor_D1 = viscosity_mat.*vel_x_par_x;
    strain_tensor_D2 = viscosity_mat.*(0.5*(vel_x_par_y+vel_y_par_x));
    strain_tensor_x = par_x(strain_tensor_D1)+par_y(strain_tensor_D2);
    strain_tensor_y = par_x(strain_tensor_D2)-par_y(strain_tensor_D1);
    
    %% Advection vector for Navier Stroke Equation
    advection_x = x_vel.*vel_x_par_x;
    advection_y = y_vel.*vel_y_par_y;
    
    %% Navier Stroke Equation:
    vel_x_par_t = ((strain_tensor_x - press_grad_x + (chem_pot_mat.*bact_frac_grad_x))./density_mat) - advection_x;
    vel_y_par_t = ((strain_tensor_y - press_grad_y + (chem_pot_mat.*bact_frac_grad_y))./density_mat) - advection_y;
    
    x_vel = x_vel + del_t*vel_x_par_t;
    y_vel = y_vel + del_t*vel_y_par_t;
    
    %% Finalizing the solvent, total bacteria volume fractions and antimicrobial concentration matrices
    for x = 0:xgridlen
        for y = 0:ygridlen
            if(bact_live_vol_frac(y+1,x+1) > 1)
                bact_live_vol_frac(y+1,x+1) = 1;
            elseif(bact_live_vol_frac(y+1,x+1)< 0)
                bact_live_vol_frac(y+1,x+1) = 0;
            end
        end
    end
    
    for x = 0:xgridlen
        for y = 0:ygridlen
            if(bact_dead_vol_frac(y+1,x+1) > 1)
                bact_dead_vol_frac(y+1,x+1) = 1;
            elseif(bact_dead_vol_frac(y+1,x+1)< 0)
                bact_dead_vol_frac(y+1,x+1) = 0;
            end
        end
    end
    
    bact_vol_frac = bact_live_vol_frac + bact_dead_vol_frac;
    for x = 0:xgridlen
        for y = 0:ygridlen
            if(bact_vol_frac(y+1,x+1) > 1)
                bact_vol_frac(y+1,x+1) = 1;
            elseif(bact_vol_frac(y+1,x+1)< 0)
                bact_vol_frac(y+1,x+1) = 0;
            end
        end
    end
    
    
    sol_vol_frac = 1 - bact_vol_frac;
    antimicrobe_con_mat = sol_con./sol_vol_frac;
    
    bact_live_vol_frac(isnan(bact_live_vol_frac))==0;
    bact_dead_vol_frac(isnan(bact_dead_vol_frac))==0;
    bact_vol_frac(isnan(1-bact_live_vol_frac))==0;
    
    mean_live = [mean_live, mean(nonzeros(bact_live_vol_frac), 'all')];
    mean_dead = [mean_dead, mean(nonzeros(bact_dead_vol_frac), 'all')];
    mean_bact = [mean_bact, mean(nonzeros(bact_vol_frac), 'all')];
    time_arr = [time_arr,(counter*del_t)];
    
    
end

mean_live(isnan(mean_live))==0;
mean_dead(isnan(mean_dead))==0;
mean_bact(isnan(mean_bact))==0;

%% PLOT INITIAL DEAD BACT VOL

% conc_fig = figure
% contour(in_l)
% colorbar;
% caxis([min max]);
% title("Initial Live Bacteria Volume fraction")

%% PLOT FINAL LIVE BACT VOL

% conc_fig = figure
% contour(bact_live_vol_frac)
% colorbar;
% caxis([min max]);
% title("Final Live Bacteria Volume fraction")

%% PLOT FINAL DEAD BACT VOL

conc_fig = figure
contour(bact_dead_vol_frac)
title("Final Dead Bacteria Volume fraction")
colorbar; caxis([min max]);

%% PLOT CONC PROFILE

% xx = 0:1:end_time;
% splive = spline(time_arr,mean_live,xx);
% spdead = spline(time_arr,mean_dead,xx);
% conc_fig = figure
% title("Mean Volume fraction for bacteria")
% % subtitle("[NisinA] = " + antimicrobe_con_d)
% subtitle("[NisinPV] = " + antimicrobe_con_d)
% ylim([0 1])
% hold on
% plot(time_arr, mean_live,'o', xx,splive, 'LineWidth', 2)
% plot(time_arr, mean_dead,'o',xx,spdead,  'LineWidth', 2)
% % plot(time_arr, mean_bact, "LineWidth", 2.5)
% % legend("live", "dead", "total")
% legend("","Live Bact Vol frac(w/o biocide)", "", "Dead Bact Vol frac (with Nisin PV)", 'Location', 'Best')
% xlabel(" Time (min)")
% ylabel("Volume fraction of bacteria")
% hold off



end

