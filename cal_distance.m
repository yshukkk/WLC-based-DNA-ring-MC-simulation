function distance = cal_distance(ads, N, theta)

x = 0;
y = 0;

for k = 1:N
    x = x + ads*cos(k*theta);
    y = y + ads*sin(k*theta);
end

distance = sqrt(x^2+y^2);

end