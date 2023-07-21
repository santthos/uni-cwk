clc
clear

increments = 1000;
frequency_max = 2 * pi * 14.4;
frequency = 0:frequency_max/increments:frequency_max;
frequency_value = 0;


for count = 1:increments+1;
    frequency_value = frequency(1,count);
    z(1,count) = (15 * frequency_value^2 * 2.9 * 10^-3)/ sqrt((18400-(15 * frequency_value^2))^2 + (580 * frequency_value)^2);
    count = count +1;
end



increments_2 = 1000;
w_max = 14.4 * 2 * pi;
w = 0:w_max/increments_2:w_max;
w_value = 0;
wn = sqrt(18400/15);

for count_2 = 1:increments_2 + 1;
    w_value = w(1,count_2);
    acceleration(1,count_2) = ((sqrt(18400^2 + (580 * w_value)))/(sqrt((18400 - (15 * w_value^2))^2 + (580 * w_value)^2)) * (w_value/wn)^2 * (2.9 * 10^-3 * wn^2));
    count_2 = count_2 +1;
end     

increments_3 = 1000;
r_max = w_max/wn;
r = 0:r_max/increments_3:r_max;
frequency_value = 0;

increments_4 = 10;


figure;
subplot 411;
plot(frequency, z)
subplot 412;
plot(w, acceleration)

for count_4 = 1:increments_4;
    
    c = count_4 * 100;
    
for count_3 = 1:increments_3 + 1;
    frequency_value = frequency(1,count_3);
    z_y(1,count_3) = ((frequency_value/wn)^2)/ sqrt(((1-((frequency_value/wn)^2))^2) + (2 * (c/(2*sqrt(15 * 18400))) * (frequency_value/wn))^2);
    count_3 = count_3 +1;
    
end

subplot 413;
plot(r, z_y)
hold on 

end
hold off

increments_5 = 1000

for count_6 = 1:increments_4;
    
    c = count_6 * 100;
    
for count_5 = 1:increments_5 +1;
    accel_Ywn(1,count_5) = ((sqrt(18400^2 + (c * w_value)))/(sqrt((18400 - (15 * w_value^2))^2 + (c * w_value)^2)) * (r(1,count_5))^2 );
    count_5 = count_5 +1;
   
end

subplot 414;
plot(r, accel_Ywn)
hold on

end

hold off





