%% Identify break and max dissipation points
function all_diagnos = find_diss_peaks

load('all_diagnos.mat');
inds = find(all_diagnos.diagnos.Time>20); % Index to remove the initial gate removal mixing from the diagnostics

Max_diss = all_diagnos.diagnos.Max_diss(inds);
smoothed_MaxDiss = smoothdata(Max_diss);
Time = all_diagnos.diagnos.Time(inds);

%%

subplot(2, 1, 1)
plot(Time, smoothed_MaxDiss)
hold on
plot(Time, Max_diss)
subplot(2, 1, 2)
plot(Time, smoothed_MaxDiss./Time.^1.1)

[~, ind] = max(smoothed_MaxDiss./Time.^2); 
MaxDiss_tot = Max_diss(ind);

%%
all_diagnos.Energetics.MaxDiss_time = round(Time(ind));
all_diagnos.Energetics.MaxDiss_tot = MaxDiss_tot;
%[MaxDiss_tot, ind] = max(Max_diss(inds));

save('all_diagnos.mat', 'all_diagnos')
end