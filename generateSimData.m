% 2 sigmas
% 2 hazard rates
% 5 simulated participant hazard rates
% 10 batches for each set of parameters
% 3 trial lengths
% 2*2*5*10*3 = 600 simulations


sigma = [30,110,160];
hTrue = [0.05,0.95];
hSubj = 0.05:.05:0.95;
tick = 0;
for index = 0:9
    for i = 1:length(sigma)
        for j = 1:length(hTrue)
            for k = 1:length(hSubj)
                tick = tick+1
                createStruct_Sim(sigma(i), hTrue(j), hSubj(k), index);
            end
        end
    end
end

subjids = cell(0,1);
m = [50,100,150];
tick = 1;
for index = 0:9
    for h = 1:length(m)
        for i = 1:length(sigma)
            for j = 1:length(hTrue)
                for k = 1:length(hSubj)
                    
                    str_m = num2str(m(h), '%03.0f');
                    str_sigma = num2str(sigma(i), '%03.0f');
                    str_hTrue = num2str(hTrue(j)*100, '%02.0f');
                    str_hSubj = num2str(hSubj(k)*100, '%02.0f');
                    str_index = num2str(index, '%02.0f');
                    
                    subjids{tick,1} = ['sim' str_m str_sigma str_hTrue str_hSubj str_index];
                    
                    tick = tick+1;
                end
            end
        end
    end
end


save('Data/Simulated_Analysis/subj_info.mat','subjids')