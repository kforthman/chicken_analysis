m = [50,100,150];
sigma = [110,160];
hTrue = [0.05,0.95];
hSubj = [0.05,0.30,0.50,0.70,0.95];
tick = 0;
for index = 0:9
    for h = 1:3
        for i = 1:2
            for j = 1:2
                for k = 1:5
                    tick = tick+1
                    createStruct_Sim(m(h), sigma(i), hTrue(j), hSubj(k), index)
                end
            end
        end
    end
end

tick = 1;
subj_info = zeros(600,1);
for index = 0:9
    for h = 1:3
        for i = 1:2
            for j = 1:2
                for k = 1:5
                    subj_info(tick,1) = cellstr(sprintf('sim_%d_%d_%.2f_%.2f_%d.mat', m(h), sigma(i), hTrue(j), hSubj(k), index))
                    tick = tick+1
                end
            end
        end
    end
end