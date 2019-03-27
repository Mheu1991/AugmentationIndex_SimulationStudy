% Simulate a batch of CircAdapt simulations from take-off model imposing
% following:
        % - maximal shortening velocity of left ventricle sarcomeres 
        % - tube stiffness      
% Maarten Heusinkveld 16-03-2019

 addpath('..\CircAdapt Model\')
 mkdir('..\','Simulations');

 % Initialisation       
 vMax_range       =  2:1:10;                % replace vMax for patches Lv1 and Sv1
 tube_stiff       = -2:1:14;                % add constant to stiffness of vascular tree

 initial_sim_time = 60;
 mmHgToPa         = 133.33;                 % 1 mmHg = 133.33 Pa
      
 try
 
     for i = 1:length(vMax_range)

         for j = 1:length(tube_stiff)

            load('PInitNEW.mat')            
           
            P.General.p0 = 92*mmHgToPa;
            
            Put('Patch','vMax',{'Lv1','Sv1'},vMax_range(i))
            Put('Patch','vMax',{'La1','Ra1'},14) % double vMax of atria so that it is equal to the CircAdapt version of Walmsley et al. 2015 PLOS CB
            
            k_init_tube       = Get('Tube','k','All');
            k_init_artven     = Get('ArtVen','k','All');
            k_updated_tube    = k_init_tube + repmat(tube_stiff(j),1,length(k_init_tube));
            k_updated_artven  = k_init_artven + repmat(tube_stiff(j),2,length(k_init_artven));
            Put('Tube','k','All',k_updated_tube)
            Put('ArtVen','k','All',k_updated_artven)

            P.General.DtSimulation  = initial_sim_time * P.General.tCycle;   
            
            save P P
            
            CircAdapt    

            P.General.DtSimulation = 4 * P.General.tCycle;
            P.General.Dt           = 1e-3;
            P.General.AdaptFunction = 'Adapt0P';        
            P.General.AdaptFeedback = 0;
            
            save P P
            
            CircAdapt
            CircDisplay

            filename = ['..\Simulations\PNew','vMax',num2str(vMax_range(i)),'kTube',num2str(tube_stiff(j)),'.mat'];
            save(filename,'P');            
            
         end       
      end
 catch
 end

