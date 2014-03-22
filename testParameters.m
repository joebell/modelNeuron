%% testParameters.m
%  
%  This shows an example of how to use the simulation to test a range of
%  parameters.  It varies the excitatory pre-synaptic rate over a range of
%  10 Hz to 40 Hz, and measures the post-synaptic rate and the CV of the
%  interspike intervals.  Then it plots them.
%
%  - JSB & AEB 3/2013
function testParameters()

    exRates = [10,15,20,25,30,35,40];   % A list rates to test (Hz)
    stepSize = .0001;                   % Simulation step size (sec)
    convergenceTime = 400;              % Allow convergence (sec)
    testTime =        100;              % Time over which to measure (sec)
    
    % For each excitatory rate in the list
    for exRateN = 1:length(exRates);
        
        % Display what rate we're testing
        disp(exRates(exRateN));
        
        % Create a model neuron
        aNeuron = modelNeuron;
        % Modify properties we care about from the defaults
        aNeuron.exSynapses.rate = exRates(exRateN);
        aNeuron.exSynapses.Aplus    =      .005; % Magnitude of synapse strengthening
        aNeuron.exSynapses.Aminus   = 1.05*.005; % Magnitude of synapse weakening
        
        % Step simulation enough that we've reached convergence
        for n = 1:round(convergenceTime/stepSize)
            aNeuron.stepTime(stepSize);
        end
        
        % Run the simulation for a time and collect a spike raster
        rasterTrace = [];
        for n = 1:round(testTime/stepSize)
            aNeuron.stepTime(stepSize);
            rasterTrace(end+1) = aNeuron.spike;
        end
        
        % Measure the mean rate and CV, store them in an array
        meanRate(exRateN) = nnz(rasterTrace)/testTime;
        spikeNs = find(rasterTrace > 0);    % Find the spike sample #'s
        spikeNDiffs = diff(spikeNs);        % # of samples between spikes
        ISIs = spikeNDiffs.*stepSize;       % Convert to ISI (sec)
        CV(exRateN) = std(ISIs)./mean(ISIs);% Calculate the CV
        
    end % End for each excitatory rate to test
    
    % Plot the mean rates and CVs for each excitatory rate tested
    subplot(1,2,1);
    plot(exRates,meanRate,'-o');
    xlabel('Pre-synaptic rate');
    ylabel('Post-synaptic rate'); 
    ylim([0 40]);
    
    subplot(1,2,2);
    plot(exRates,CV,'-o');
    xlabel('Pre-synaptic rate');
    ylabel('CV of ISI'); 
    ylim([0 1]);
    
        
        
        
