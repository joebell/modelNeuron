%%      modelNeuron.m
%
%  This class defines a model neuron with STDP, as described in 
%  Song, et. al., 2000.  We find it should be run with time steps of 
%  .1 ms or smaller, in order to reproduce the key findings.
%
%
%  Properties of the simulation are defined in the constructor.
%
%
%  In our hands, depending on the input, the simulation takes 
%  at least 100 seconds of simulated time to converge on a stable 
%  post-synaptic rate and distribution of conductances, and sometimes 
%  substantially longer.  Tip: Increasing the strength of STDP and the 
%  presynaptic rates will speed convergence.
%
%  To run the simulation, first construct an instance of this class:
%
%       myModelNeuron = modelNeuron();
%
%  Simulation time can then be advanced by invoking the stepTime() method, 
%  which takes the size of the desired time step (in sec.) as an argument:
%
%       myModelNeuron.stepTime(.0001);
%
%  You can access properties in a similar manner.  For example:
%
%       currentVm = myModelNeuron.Vm;
%       isSpiking = myModelNeuron.spike;
%
%
%
% - JSB and AEB, 2/2013
classdef modelNeuron < handle
   
   %% Properties Block: 
   %   
   %  Object properties are declared here, but default values are defined
   %  below, in the constructor.
   properties 
      
       Vm     
       gEx    
       gIn    
       spike   
              
       tauM   
       Vrest  
       Eex     
       Ein    
       tauEx   
       tauIn      
       Vthresh
       Vreset 
       
       Nex
       Nin
       exSynapses
       inSynapses 
       
   end % end properties
   
   
   
   
   
   %% Methods Block:
   methods
       
       
       %% Constructor
       %
       %  This creates the default model neuron.  The 'this' syntax looks 
       %  goofy, but it lets MATLAB know that we're defining the parameter
       %  values for the particular modelNeuron that we create when we 
       %  invoke the constructor.
       %
       function this = modelNeuron()

           % Variables - these change over the simulation       
           this.Vm      =   -58;    % Membrane voltage, mV   
           this.gEx     =     0;    % Total excitatory conductance
           this.gIn     =     0;    % Total inhibitory conductance
           this.spike   = false;    % 'true' if the neuron is spiking, else 'false'

           % Neuron properties
           this.tauM    =  .020;    % Membrane time constant, sec
           this.Vrest   =   -70;    % Resting membrane voltage, mV
           this.Eex     =     0;    % Excitatory reversal potential, mV
           this.Ein     =   -70;    % Inhibitory reversal potential, mV
           this.tauEx   =  .005;    % Time constant of excitatory conductances
           this.tauIn   =  .005;    % Time constant of inhibitory conductances     
           this.Vthresh =   -54;    % Spike threshold voltage, mV
           this.Vreset  =   -60;    % Post-spike reset voltage, mV   
           
           %% Define excitatory synapses                    
           this.Nex                 = 1000; % # of excitatory synapses           
           this.exSynapses.rate     =   25; % Rate of pre-synaptic APs           
           % Maximum peak excitatory conductance
           this.exSynapses.gMax     = .015;
           % Peak excitatory conductance for each synapse a, these are
           % initially set to gMax
           this.exSynapses.gA       = this.exSynapses.gMax.*ones(this.Nex,1);
           % A vector listing which presynaptic synapses are firing
           this.exSynapses.preSynapticSpike = zeros(this.Nex,1);
           % M and Pa are house keeping variables for implementing STDP
           % rule.  Both decay exponentially toward zero.
           %    Pa is positive, used to increase the strength of 
           %       synapses, and is incremented on pre-synaptic spiking.
           %       On postsynaptic spiking, gA -> gA + P*gMax
           %    M is negative, used to decrease the strength of synapses, 
           %       and is decremented on post-synaptic spiking.
           %       On presynaptic spiking,  gA -> gA + M*gMax          
           this.exSynapses.Pa       = zeros(this.Nex,1);
           this.exSynapses.M        = 0;          
           % STDP parameters for excitatory synapses
           this.exSynapses.Aplus    =      .005; % Magnitude of synapse strengthening
           this.exSynapses.Aminus   = 1.05*.005; % Magnitude of synapse weakening
           this.exSynapses.tauPlus  =      .020; % Time constant of strengthening (sec)
           this.exSynapses.tauMinus =      .020; % Time constant of weakening (sec)
                     
           
           %% Define inhibitory synapses (nb. These are not subject to STDP)
           this.Nin                 = 200;  % # of inhibitory synapses             
           this.inSynapses.rate     =  10;  % Rate of presynaptic APs
           % A vector listing which presynaptic synapses are firing
           this.inSynapses.preSynapticSpike = zeros(this.Nin,1);
            % Peak excitatory conductance for each synapse a
           this.inSynapses.gA               = .050*ones(this.Nin,1);  
           
           
       end % End modelNeuron()
       
       
       
       
       
       
       %% Advance time by one step of length dT
       function stepTime(this, dT)
           
           % Determine if the neuron will spike, or not
           if (this.Vm > this.Vthresh)                  
               this.spike = true;        % Note there's a spike
               % Reset Vm to the reset voltage
               this.Vm = this.Vreset;

               % Update learning rule M
               this.exSynapses.M = this.exSynapses.M - this.exSynapses.Aminus;
               % Update conductances as a result of the learning rule applied
               % to post-synaptic spikes.
               this.exSynapses.gA = this.exSynapses.gA + this.exSynapses.Pa.*this.exSynapses.gMax;           
               % Don't allow conductances out of the range [0,gMax]
               this.exSynapses.gA = this.exSynapses.gA - ...
                    (this.exSynapses.gA > this.exSynapses.gMax).*...
                    (this.exSynapses.gA - this.exSynapses.gMax);
               this.exSynapses.gA = this.exSynapses.gA - ...
                    (this.exSynapses.gA < 0).*...
                    (this.exSynapses.gA);           
           else % If it doesn't spike...          
               this.spike = false;       % Note there's no spike
               % Update membrane voltage based on conductances
               dV = (dT/this.tauM)*(this.Vrest - this.Vm ...
                                    + this.gEx*(this.Eex - this.Vm) ...
                                    + this.gIn*(this.Ein - this.Vm) );
               this.Vm = this.Vm + dV;
           end
           
           % Allow conductances to decay exponentially
           dgEx = -this.gEx*dT/this.tauEx;
           dgIn = -this.gIn*dT/this.tauIn;
           this.gEx = this.gEx + dgEx;
           this.gIn = this.gIn + dgIn;
           
           % Generate Poisson presynaptic spikes, 1 for spike, 0 for none
           this.exSynapses.preSynapticSpike = (rand(this.Nex,1) < dT*this.exSynapses.rate);
           this.inSynapses.preSynapticSpike = (rand(this.Nin,1) < dT*this.inSynapses.rate);
           
           % Presynaptic spikes generate conductances in the post-synaptic
           % cell
           exCond = this.exSynapses.preSynapticSpike.*this.exSynapses.gA;
           inCond = this.inSynapses.preSynapticSpike.*this.inSynapses.gA;
           this.gEx = this.gEx + sum(exCond);         
           this.gIn = this.gIn + sum(inCond);
           
           % Update learning rule: Pa increases conductances on
           % post-synaptic spiking, and is incremented on pre-synaptic
           % spiking.
           this.exSynapses.Pa = this.exSynapses.Pa +...
                   this.exSynapses.preSynapticSpike.*this.exSynapses.Aplus;
           
           % Update the conductances as a result of the learning rule
           % applied to pre-synaptic spikes.
           this.exSynapses.gA = this.exSynapses.gA + ...
               this.exSynapses.preSynapticSpike.*this.exSynapses.M*this.exSynapses.gMax;           
           % Don't allow conductances out of the range [0,gMax]
           this.exSynapses.gA = this.exSynapses.gA - ...
                (this.exSynapses.gA > this.exSynapses.gMax).*...
                (this.exSynapses.gA - this.exSynapses.gMax);
           this.exSynapses.gA = this.exSynapses.gA - ...
                (this.exSynapses.gA < 0).*...
                (this.exSynapses.gA);         
                                
           % The learning rule functions M and Pa decay exponentially
           dM = -this.exSynapses.M*dT/this.exSynapses.tauMinus;
           this.exSynapses.M = this.exSynapses.M + dM;
           dPa = -this.exSynapses.Pa.*dT/this.exSynapses.tauPlus;
           this.exSynapses.Pa = this.exSynapses.Pa + dPa;
                   
       end % End stepTime()
       
       
       
       
       
       
   end % End methods
end % End classdef