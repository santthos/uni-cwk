% the data below is the raw data collected along with the relevant unit
% conversions
Fq = [1.2566 1.8850 2.5133 3.1416 3.7699 4.3982];

Ar = [0.966 0.937 0.942 0.923 0.969 0.923];

Ph = [-0.1074 -0.2096 -0.3142 -0.4188 -0.2611 -0.4438];

Ph2 = [-6.156 -12.0108 -18 -23.994 -14.9605 -25.42909];

Ar2 = [-0.3005 -0.5652 -0.5190 -0.6960 -0.2735 -0.6960];

% this part of the code produces the bode plot estimation 
zfr = Ar.*exp(1i*Ph);
data = frd(zfr,Fq);
% by altering the nz and np values a transfer function which fit the
% experimental data was found
nz = 5; np = 6;
Gc = tfest(data,np,nz);

% this part of the code takes the amplitude ratios from the bode plot and
% converts it into dB units
[mag, phase, wout] = bode(Gc);
mag2 = 20 * log10(mag(:,:));

% below is the part of the code which plots the experimental data onto the
% bode plot
subplot(2,1,1)
figure1 = semilogx(wout,mag2);
hold on
plot(Fq,Ar2,'x');
title('Amplitude Ratio vs Frequency');
xlabel('Frequency (Rad/s)') 
ylabel('Amplitude Ratio dB') 
hold off

phase2 = phase(:,:);

subplot(2,1,2)
figure2 = semilogx(wout,phase2);
title('Phase Angle vs Frequency')
hold on
plot(Fq,Ph2,'x');
xlabel('Frequency (Rad/s)') 
ylabel('Phase Angle (Degrees)') 
hold off

% this part simply shows the transfer function which fits with the
% experimental data
Gc


