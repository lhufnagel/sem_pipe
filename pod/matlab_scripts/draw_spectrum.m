if ~exist('time_coeffs')
  load('time_coeffs_asym.mat');
end

figure(1)
clf;
hold on;

modes = [2,3,4,5];


for mode = modes
  [pxx,f] = pwelch(time_coeffs(mode,:),150,75,1000,4);

  %semilogy(f,pxx)
  if (mode > 3)
   plot(f,10*log10(pxx),'x-')
  else
    plot(f,10*log10(pxx))
  end

  xlabel('Strouhal number')
  ylabel('Magnitude (dB)')
end
grid on


legend({'Time coef. 1', 'Time coef. 2','Time coef. 3','Time coef. 4'},'Location','SOUTHWEST')
set(gca,'fontsize', 18);
%legend(num2str(modes'-1))

% create smaller axes in top right, and plot on it
axes('Position',[.5 .5 .4 .42])
box on
hold on

for mode = modes
  [pxx,f] = pwelch(time_coeffs(mode,:),500,250,2500,4);

  %semilogy(f,pxx)
  if (mode > 3)
   plot(f,10*log10(pxx),'x-')
  else
    plot(f,10*log10(pxx))
  end

  xlabel('Strouhal number')
  ylabel('Magnitude (dB)')
end
grid on
axis([0.05 0.4 -50 -10])

set(gca,'fontsize', 18);
