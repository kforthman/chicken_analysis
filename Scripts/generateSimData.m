% 2 sigmas
% 2 hazard rates
% 19 simulated participant hazard rates
% 6 lapse rates
% 13 noise rates
% 10 batches for each set of parameters
% 3 trial lengths
% 2*2*19*6*13*10*3 = 177,840 simulations


sigma = [110,160];
hTrue = [0.05,0.95];
hSubj = 0.05:.05:0.95;
lapse = 0:0.05:.25;
noise = [0.01,0.5:0.5:6];
tick = 0;
for index = 0:9
    for i = 1:length(sigma)
        for j = 1:length(hTrue)
            for k = 1:length(hSubj)
                for m = 1:length(lapse)
                    for n = 1:length(noise)
                        createStruct_Sim(sigma(i), hTrue(j), hSubj(k), lapse(m), noise(n), index);
                        tick = tick+1
                        tick = tick+1
                        tick = tick+1
                    end
                end
            end
        end
    end
end

subjids = cell(0,1);
N = [50,100,150];
tick = 1;

for h = 1:length(N)
    for i = 1:length(sigma)
        for j = 1:length(hTrue)
            for k = 1:length(hSubj)
                for m = 1:length(lapse)
                    for n = 1:length(noise)
                        for index = 0:9
                            
                            str_m = num2str(N(h), '%03.0f');
                            str_sigma = num2str(sigma(i), '%03.0f');
                            str_hTrue = num2str(hTrue(j)*100, '%02.0f');
                            str_hSubj = num2str(hSubj(k)*100, '%02.0f');
                            str_lapse = num2str(lapse(m)*100, '%02.0f');
                            str_noise = num2str(noise(n)*100, '%03.0f');
                            str_index = num2str(index, '%02.0f');
                            
                            subjids{tick,1} = [...
                                'sim-' ...
                                'N_'      str_m ...
                                '-sigma_' str_sigma ...
                                '-hT_'    str_hTrue ...
                                '-hS_'    str_hSubj ...
                                '-lapse_' str_lapse ...
                                '-noise_' str_noise ...
                                '-' str_index ...
                                ];
                            
                            tick = tick+1;
                            
                        end
                    end
                end
            end
        end
    end
end


save('../Output/Simulated_nl_Analysis/subj_info.mat','subjids')