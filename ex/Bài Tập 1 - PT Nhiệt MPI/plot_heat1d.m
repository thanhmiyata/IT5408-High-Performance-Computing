data = readmatrix("T_serial.csv"); % hoặc T_omp.csv / T_mpi.csv
t = data(:,1);
T = data(:,2:end);

M = size(T,2);
x = linspace(0,1,M);

figure; hold on;
idx = round(linspace(1,length(t),5)); % 5 thời điểm
for k = 1:length(idx)
    plot(x, T(idx(k),:));
end
xlabel('x');
ylabel('Temperature (C)');
legend(string(t(idx)));
grid on;