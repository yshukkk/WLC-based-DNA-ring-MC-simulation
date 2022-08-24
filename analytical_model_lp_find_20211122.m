N = 30;
ring = 48:100;
ring(length(ring)+1) = 1000000;
ring = [48, 56, 64, 70, 80, 100, 10000]; % used ring ssDNA size and linear DNA

exp_ring_30 = [48, 56, 64, 70, 80, 100]; % add data of dye-to-dye distance
exp_30 = [0.77202, 0.87058, 0.92873, 0.95364, 0.98067, 0.99527];
exp_30_dR = [-1.841427987, -1.013385007, -0.524854276, -0.315598739, -0.088498535, 0.034159966];
exp_30_dR_sigma = [0.054552189, 0.048986291, 0.057464159, 0.046836324, 0.061206428, 0.063613845];

exp_ring_32 = [48, 56, 64, 70, 80, 100];
exp_32_dR = [-1.725243855, -0.765576093, -0.444043213, -0.261878677, -0.139838498, -0.009510085];
exp_32 = [0.78787, 0.84862, 0.90222, 0.94891, 0.96806, 0.98028, 0.98844];
exp_32_dR_sigma = [0.01300516, 0.019335654, 0.020243874, 0.012150833, 0.007535739, 0.025133264];

exp_ring_34 = [48, 56, 64, 70, 80, 100];
exp_34_dR = [-1.757865224, -0.732851628, -0.416707315, -0.225867311, -0.021501952, 0.133841458];
exp_34_dR_sigma = [0.054015237, 0.070790672, 0.053660088, 0.05961084, 0.086806022, 0.084120699];

exp_ring_35 = [48, 56, 64, 70, 80, 100];
exp_35_dR = [-2.01456, -0.88471, -0.53629, -0.34803, -0.17649, -0.0473];
exp_35_dR_sigma = [0.02744, 0.02254, 0.02712, 0.0253, 0.02565, 0.05807];

exp_ring_36 = [48, 56, 64, 70, 80, 100];
exp_36_dR = [-1.763303762, -0.63924159, -0.413275313, -0.329122342, -0.125489139, -0.07888166];
exp_36_dR_sigma = [0.042670177, 0.033181958, 0.040647872, 0.031689472, 0.034160386, 0.032165669];

exp_ring_38 = [56, 64, 70, 80, 100];
exp_38_dR = [-1.198877248, -0.211244765, -0.216362238, -0.119504747, -0.02766357];
exp_38 = [0.96185, 0.96123, 0.97302, 0.9842];
exp_38_dR_sigma = [0.117122755, 0.04828268, 0.051669804, 0.06976874, 0.048326731];

exp_ring_40 = [56, 64, 70, 80, 100];
exp_40_dR = [-1.849174719, -0.371771763, -0.263753937, -0.003072965, 0.005409477];
exp_40_dR_sigma = [0.239783627, 0.029453454, 0.031878952, 0.043563528, 0.043929931];

n_0 = 18;
theta = 0;
ads = 0.34;
ass = 0.7;
Lss = n_0*ass;
Lpds = 15;
Lpds_flex = 10;
Lpss = 3;
bending_module = Lpds/ads;
bending_module_flex = Lpds_flex/ads;
bending_module_end = 2*Lpds_flex*Lpss/(ads*Lpss+ass*Lpds_flex);
resolution = 1000;
a = 1;

exp_x = exp_ring_35;
exp_y = exp_35_dR;
exp_sigma = exp_35_dR_sigma;

ss_size = exp_x-N;
%ss_size(length(ss_size)+1) = 1000000;
interval = 0.1;

energy_bend = zeros(1,resolution);
energy_end_bending = zeros(1,resolution);
y=0;

for j = 1:interval:100
    Lpds = j;
    bending_module = Lpds/ads;
    y = y+1;
    aver_dye = 0;
    diff_aver_dye = 0;
    
    for i = ss_size
        n_0 = i;
        Lss = n_0*ass;
        
        for k = 1:resolution
            theta(k) = 2*pi/(N)*k/resolution;
            %energy_bend(k) = (N-1)*(bending_module*(1-cos(theta(k))))^a; % cosine potential
            %energy_bend(k) = 29*bending_module/2*theta(k); %SLEC potential
            energy_bend(k) = bending_module/2*(N)*theta(k)^2; %harmonic potential
            %energy_bend(k) = bending_module/2*(N-11)*theta(k)^2 + bending_module_flex/2*10*theta(k)^2;
            
            distance(k) = cal_distance(ads,N,theta(k));
            distance_dye(k) = cal_distance(ads,24,theta(k));
            %energy_end_bending(k) = 2*bending_module_end*(1-end_angle_cos(ads,N,theta(k)));
            inver_lang(k) = distance(k)/Lss*(3-2.7*distance(k)/Lss+0.7*(distance(k)/Lss)^2)/((1-distance(k)/Lss)*(1+0.1*distance(k)/Lss));
            
            if 0.95*Lss <= distance(k)
                energy_stretching(k) = inf;
            else
                energy_stretching(k) = n_0*(distance(k)/Lss*inver_lang(k)+log(inver_lang(k)/sinh(inver_lang(k))));
                %energy_stretching(k) = n_0*(-7/2*(distance(k)/Lss)^2 + 89*distance(k)/Lss - log(1-distance(k)/Lss) - 900*log(10+distance(k)/Lss) + 900*log(10));
                force_stretching(k) = n_0/Lss*inver_lang(k);
            end
        end
        
        fret = 1./(1+(distance_dye/6.5).^6);
        diff_dye = distance_dye-max(distance_dye);
        
        total_energy = energy_bend + energy_stretching;
        total_energy(imag(total_energy) ~= 0) = 10000;
        
        
        [energy,angle_k] = min(total_energy);
        cal_dye(find(ss_size == i)) = distance_dye(angle_k);
        cal_diff_dye = cal_dye-max(cal_dye);
        cal_fss(find(ss_size == i)) = force_stretching(angle_k);
        cal_fss(find(ss_size == i)) = n_0/Lss*distance_dye(angle_k)/Lss*(3-2.7*distance_dye(angle_k)/Lss+0.7*(distance_dye(angle_k)/Lss)^2)/((1-distance_dye(angle_k)/Lss)*(1+0.1*distance_dye(angle_k)/Lss));
    end
    least_sq(y) = sum((exp_y - cal_diff_dye).^2./exp_sigma.^2); %weight sigma
    %least_sq(y) = sum((exp_y - diff_aver_dye).^2);
    func(y,:) = cal_diff_dye;
    if y >= 3
        dfunc_per_dLp(y-1,:) = (func(y,:)-func(y-2,:))./(2*interval*exp_sigma);
    end
    
    disp(j);
end

[min_least,lp_fit] = min(least_sq);

beta = dfunc_per_dLp(lp_fit,:);
error = sqrt(min_least/(length(ss_size)-2)/sum(beta.^2))
lp_sol = lp_fit*interval + 0.9
plot(1:0.1:100,least_sq);

bending_module = 50/ads;

for i = 1:100
    n_0 = i;
    Lss = n_0*ass;
    
    for k = 1:resolution
        theta(k) = 2*pi/(N)*k/resolution;
        %energy_bend(k) = (N-1)*(bending_module*(1-cos(theta(k))))^a; % cosine potential
        %energy_bend(k) = 29*bending_module/2*theta(k); %SLEC potential
        energy_bend(k) = bending_module/2*(N)*theta(k)^2; %harmonic potential
        %energy_bend(k) = bending_module/2*(N-11)*theta(k)^2 + bending_module_flex/2*10*theta(k)^2;
        
        distance(k) = cal_distance(ads,N,theta(k));
        distance_dye(k) = cal_distance(ads,24,theta(k));
        %energy_end_bending(k) = 2*bending_module_end*(1-end_angle_cos(ads,N,theta(k)));
        inver_lang(k) = distance(k)/Lss*(3-2.7*distance(k)/Lss+0.7*(distance(k)/Lss)^2)/((1-distance(k)/Lss)*(1+0.1*distance(k)/Lss));
        
        if 0.95*Lss <= distance(k)
            energy_stretching(k) = inf;
        else
            energy_stretching(k) = n_0*(distance(k)/Lss*inver_lang(k)+log(inver_lang(k)/sinh(inver_lang(k))));
            %energy_stretching(k) = n_0*(-7/2*(distance(k)/Lss)^2 + 89*distance(k)/Lss - log(1-distance(k)/Lss) - 900*log(10+distance(k)/Lss) + 900*log(10));
            force_stretching(k) = n_0/Lss*inver_lang(k);
        end
    end
    
    fret = 1./(1+(distance_dye/6.5).^6);
    diff_dye = distance_dye-max(distance_dye);
    
    total_energy = energy_bend + energy_stretching;
    total_energy(imag(total_energy) ~= 0) = 10000;
    
    
    [energy,angle_k] = min(total_energy);
    cal_dye(i) = distance_dye(angle_k);
    cal_diff_dye = cal_dye-max(cal_dye);
    cal_fss(i) = force_stretching(angle_k);
    cal_fss2(i) = n_0/Lss*distance_dye(angle_k)/Lss*(3-2.7*distance_dye(angle_k)/Lss+0.7*(distance_dye(angle_k)/Lss)^2)/((1-distance_dye(angle_k)/Lss)*(1+0.1*distance_dye(angle_k)/Lss));
end
plot([1:100]+N, cal_diff_dye, '--', exp_x, exp_y, 'o');
%axis([40 100 -3 0])
