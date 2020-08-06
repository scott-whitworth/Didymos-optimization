%This function is designed to take Power, efficiency, uExhaust, and mdot
%and put them all into one plot. For a reference on what we are trying to
%produce look at the 2010 Sankaran paper, Fig.1. We have taken data off of
%three other papers to give us better data for each thruster.


function [ ] = EtaPowerGraphs( thrusterID )

INCREMENT_CONSTANT = 20;
%%%%%%%%%%%%%%%% BHT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if thrusterID == 1
    
      
    thrusterName = 'BHT';
    powerInPlotting = xlsread('ThrustersData.xlsx','BHT','b33:b114');
    etaPlotting = xlsread('ThrustersData.xlsx','BHT','G33:G114');
  
    %first set with mDot = 5E-05 kg/s
    powerIn1 = xlsread('ThrustersData.xlsx','BHT', 'A3:A6');
    eta1 = xlsread('ThrustersData.xlsx','BHT', 'F3:F6');
    
    
    %second set with mDot = 4.5E-05 kg/s
    powerIn2 = xlsread('ThrustersData.xlsx','BHT', 'A7:A10');
    eta2 = xlsread('ThrustersData.xlsx','BHT', 'F7:F10');
   
    
   %third set with mDot = 4E-05 kg/s
    powerIn3 = xlsread('ThrustersData.xlsx','BHT', 'A11:A15');
    eta3 = xlsread('ThrustersData.xlsx','BHT', 'F11:F15');
  
     
    %4th set with mDot = 3.5E-05 kg/s
    powerIn4 = xlsread('ThrustersData.xlsx','BHT', 'A16:A19');
    eta4 = xlsread('ThrustersData.xlsx','BHT', 'F16:F19');
    
    
    %fith set with mDot = 3E-05 kg/s
    powerIn5 = xlsread('ThrustersData.xlsx','BHT', 'A20:A23');
    eta5 = xlsread('ThrustersData.xlsx','BHT', 'F20:F23');
    
    
    %6th set with mDot = 2.5E-05 kg/s
    powerIn6 = xlsread('ThrustersData.xlsx','BHT', 'A24:A25');
    eta6 = xlsread('ThrustersData.xlsx','BHT', 'F24:F25');
  
    
    %7th set with mDot = 2E-05 kg/s
    powerIn7 = xlsread('ThrustersData.xlsx','BHT', 'A26:A27');
    eta7 = xlsread('ThrustersData.xlsx','BHT', 'F26:F27');
   
    
    %8th set with mDot = 1.5E-05 kg/s
    powerIn8 = xlsread('ThrustersData.xlsx','BHT', 'A28:A28');
    eta8 = xlsread('ThrustersData.xlsx','BHT', 'F28:F28');
    
       
    %%%%%%%%%%%%%%%%%%%%%%%%% ETA'S %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %this calculates the eta of the constant mdot
    
    
    
    
    %This is the line fit for when the mdot is constant
    
    pIn1 =  powerIn1 - min(powerIn1);
     pIn2 =  powerIn2 - min(powerIn2);
      pIn3 =  powerIn3 - min(powerIn3);
       pIn4 =  powerIn4 - min(powerIn4);
        pIn5 =  powerIn5 - min(powerIn5);
         pIn6 =  powerIn6 - min(powerIn6);
          pIn7 =  powerIn7 - min(powerIn7);
           pIn8 =  powerIn8 - min(powerIn8);
     
    %%%%%%makeing the curve smooth by making powerSteps with more values
    stepSize = [(max(pIn1) - min(pIn1))/ INCREMENT_CONSTANT;(max(pIn2) - min(pIn2))/ INCREMENT_CONSTANT;(max(pIn3) - min(pIn3))/ INCREMENT_CONSTANT;
        (max(pIn4) - min(pIn4))/ INCREMENT_CONSTANT;(max(pIn5) - min(pIn5))/ INCREMENT_CONSTANT;(max(pIn6) - min(pIn6))/ INCREMENT_CONSTANT;
        (max(pIn7) - min(pIn7))/ INCREMENT_CONSTANT;(max(pIn8) - min(pIn8))/ INCREMENT_CONSTANT;];
    
    
    powerSteps1 = max(pIn1):-stepSize(1,1):min(pIn1);
     powerSteps2 = max(pIn2):-stepSize(2,1):min(pIn2);
      powerSteps3 = max(pIn3):-stepSize(3,1):min(pIn3);
       powerSteps4 = max(pIn4):-stepSize(4,1):min(pIn4);
        powerSteps5 = max(pIn5):-stepSize(5,1):min(pIn5);
         powerSteps6 = max(pIn6):-stepSize(6,1):min(pIn6);
          powerSteps7 = max(pIn7):-stepSize(7,1):min(pIn7);
           powerSteps8 = max(pIn8):-stepSize(8,1):min(pIn8);
               
           powerSteps1 = powerSteps1';
            powerSteps2 = powerSteps2';
             powerSteps3 = powerSteps3';
              powerSteps4 = powerSteps4';
               powerSteps5 = powerSteps5';
                powerSteps6 = powerSteps6';
                 powerSteps7 = powerSteps7';
                  powerSteps8 = powerSteps8';
    
           
    
    lineFit = polyfit(sqrt(pIn1),eta1,1);
    YMDot1 = polyval(lineFit,sqrt(powerSteps1));
       
    lineFit = polyfit(sqrt(pIn2),eta2,1);
    YMDot2 = polyval(lineFit,sqrt(powerSteps2));
       
    lineFit = polyfit(sqrt(pIn3),eta3,1);
    YMDot3 = polyval(lineFit,sqrt(powerSteps3));
    
    lineFit = polyfit(sqrt(pIn4),eta4,1);
    YMDot4 = polyval(lineFit,sqrt(powerSteps4));
    
    lineFit = polyfit(sqrt(pIn5),eta5,1);
    YMDot5 = polyval(lineFit,sqrt(powerSteps5));
    
    lineFit = polyfit(sqrt(pIn6),eta6,1);
    YMDot6 = polyval(lineFit,sqrt(powerSteps6));
    
    lineFit = polyfit(sqrt(pIn7),eta7,1);
    YMDot7 = polyval(lineFit,sqrt(powerSteps7));
    
    lineFit = polyfit(sqrt(pIn8),eta8,1);
    YMDot8 = polyval(lineFit,sqrt(powerSteps8));
    
    %%%%%set powerSteps back to measued min and max values for ploting
    
    powerSteps1 = max(powerIn1):-stepSize(1,1):min(powerIn1);
     powerSteps2 = max(powerIn2):-stepSize(2,1):min(powerIn2);
      powerSteps3 = max(powerIn3):-stepSize(3,1):min(powerIn3);
       powerSteps4 = max(powerIn4):-stepSize(4,1):min(powerIn4);
        powerSteps5 = max(powerIn5):-stepSize(5,1):min(powerIn5);
         powerSteps6 = max(powerIn6):-stepSize(6,1):min(powerIn6);
          powerSteps7 = max(powerIn7):-stepSize(7,1):min(powerIn7);
           powerSteps8 = max(powerIn8):-stepSize(8,1):min(powerIn8);
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%    KEEPING UEX CONSTANT %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
      
    %First set at Uex = 25000 m/s
    powerInB1 = xlsread('ThrustersData.xlsx','BHT', 'Q3:Q4');
    etaB1 = xlsread('ThrustersData.xlsx','BHT', 'P3:P4');
    
    %Second set at Uex = 22500 m/s
    powerInB2 = xlsread('ThrustersData.xlsx','BHT', 'Q5:Q10');
    etaB2 = xlsread('ThrustersData.xlsx','BHT', 'P5:P10');
    
    %Third set at Uex = 20000 m/s
    powerInB3 = xlsread('ThrustersData.xlsx','BHT', 'Q11:Q16');
    etaB3 = xlsread('ThrustersData.xlsx','BHT', 'P11:P16');
    
    %Fourth set at Uex = 17500 m/s
    powerInB4 = xlsread('ThrustersData.xlsx','BHT', 'Q17:Q21');
    etaB4 = xlsread('ThrustersData.xlsx','BHT', 'P17:P21');
    
    
    
    pInB1 = powerInB1-min(powerInB1);
     pInB2 = powerInB2-min(powerInB2);
      pInB3 = powerInB3-min(powerInB3);
       pInB4 = powerInB4-min(powerInB4);
    
    %make a smoother line by using more Pin values
    
    stepSizeB = [(max(pInB1) - min(pInB1))/ INCREMENT_CONSTANT;(max(pInB2) - min(pInB2))/ INCREMENT_CONSTANT;
        (max(pInB3) - min(pInB3))/ INCREMENT_CONSTANT;(max(pInB4) - min(pInB4))/ INCREMENT_CONSTANT];
    
    %%%sets powerSteps to have the min start at zero for giving it a better
    %%%line fit
    
    powerStepsB1 = min(pInB1):stepSizeB(1,1):max(pInB1);
     powerStepsB2 = min(pInB2):stepSizeB(2,1):max(pInB2);
      powerStepsB3 = min(pInB3):stepSizeB(3,1):max(pInB3);
       powerStepsB4 = min(pInB4):stepSizeB(4,1):max(pInB4);
     
     
    lineFit = polyfit(sqrt(pInB1),etaB1,1);
    
     %%% replace negative slopes w/ zero
    if lineFit(1,1) < 0
        lineFit(1,1) = 0;
    end
    
    YUex1 = polyval(lineFit,sqrt(powerStepsB1));
    
    lineFit = polyfit(sqrt(pInB2),etaB2,1);
    
     %%% replace negative slopes w/ zero
    if lineFit(1,1) < 0
        lineFit(1,1) = 0;
    end
    
    YUex2 = polyval(lineFit,sqrt(powerStepsB2));
    
    lineFit = polyfit(sqrt(pInB3),etaB3,1);
    
     %%% replace negative slopes w/ zero
    if lineFit(1,1) < 0
        lineFit(1,1) = 0;
    end
    
    YUex3 = polyval(lineFit,sqrt(powerStepsB3));
   
    lineFit = polyfit(sqrt(pInB4),etaB4,1);
    
    %%% replace negative slopes w/ zero
    if lineFit(1,1) < 0
        lineFit(1,1) = 0;
    end
    
    YUex4 = polyval(lineFit,sqrt(powerStepsB4));

    %%%%% set powerSteps back to measued min and max values for ploting
    
    
     powerStepsB1 = min(powerInB1):stepSizeB(1,1):max(powerInB1);
      powerStepsB2 = min(powerInB2):stepSizeB(2,1):max(powerInB2);
       powerStepsB3 = min(powerInB3):stepSizeB(3,1):max(powerInB3);
        powerStepsB4 = min(powerInB4):stepSizeB(4,1):max(powerInB4);
         
        
    plot(powerInPlotting, etaPlotting,'+k',powerSteps1,YMDot1,'-c',powerStepsB1, YUex1, '-.r', powerSteps2,YMDot2,'-c',powerStepsB2, YUex2, '-.r',powerSteps3,YMDot3,'-c',powerStepsB3, YUex3, '-.r',powerSteps4,YMDot4,'-c',powerStepsB4, YUex4, '-.r',powerSteps5,YMDot5,'-c',powerSteps6,YMDot6,'-c',powerSteps7,YMDot7,'-c',powerSteps8,YMDot8,'-c'); 
    
    %Labels the mDot Lines
    mDotLabels = {'50 mg/s','45 mg/s','40 mg/s','35 mg/s','30 mg/s','25 mg/s','20 mg/s','15 mg/s'};
    text(powerIn1(1,1),YMDot1(1,1),mDotLabels(1,1),'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(powerIn2(1,1),YMDot2(1,1),mDotLabels(1,2),'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(powerIn3(1,1),YMDot3(1,1),mDotLabels(1,3),'VerticalAlignment','top','HorizontalAlignment','left')
    text(powerIn4(1,1),YMDot4(1,1),mDotLabels(1,4),'VerticalAlignment','bottom','HorizontalAlignment','right')
    text(powerIn5(1,1),YMDot5(1,1),mDotLabels(1,5),'VerticalAlignment','top','HorizontalAlignment','right')
    text(powerIn6(1,1),YMDot6(1,1),mDotLabels(1,6),'VerticalAlignment','top','HorizontalAlignment','right')
    text(powerIn7(1,1),YMDot7(1,1),mDotLabels(1,7),'VerticalAlignment','top','HorizontalAlignment','right')
  %  text(powerIn8(1,1),YMDot8(1,1),mDotLabels(1,8),'VerticalAlignment','bottom','HorizontalAlignment','right')
    
    %Labels the uExhaust Lines
    uExhaustLabels = {'25000 m/s','22500 m/s','20000 m/s','17500 m/s'};
    text(powerInB1(1,1),YUex1(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,1),'VerticalAlignment','bottom','HorizontalAlignment','right')
    text(powerInB2(1,1),YUex2(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,2),'VerticalAlignment','bottom','HorizontalAlignment','right')
    text(powerInB3(1,1),YUex3(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,3),'VerticalAlignment','bottom','HorizontalAlignment','right')
    text(powerInB4(1,1),YUex4(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,4),'VerticalAlignment','bottom','HorizontalAlignment','right')
    
    title(thrusterName)
    xlabel('Power (W)')
    ylabel('Efficiency')
    axis([1000,24000,.3,.8])
    legend('Measured Values','Constant mDot','Constant uExhaust','Location','southeast')



%%%%%%%%%%%%%%%%%%%% HiPEP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif thrusterID == 2
    
    
    thrusterName = 'HiPEP';
    %measured points to plot
    powerInPlotting = xlsread('ThrustersData.xlsx','HiPEP','G3:G37');
    etaPlotting = xlsread('ThrustersData.xlsx','HiPEP','D3:D37');
        
    
    %First set of data with mDot = 4.0E-06 kg/s
    powerIn1 = xlsread('ThrustersData.xlsx','HiPEP','G3:G16');
    eta1 = xlsread('ThrustersData.xlsx','HiPEP','D3:D16');
       
    %Second set of data with mDot = 6.6E-06 kg/s
    powerIn2 = xlsread('ThrustersData.xlsx','HiPEP','G17:G27');
    eta2 = xlsread('ThrustersData.xlsx','HiPEP','D17:D27');
      
    %Third set of data with mDot = 7.0E-06 kg/s
    powerIn3 = xlsread('ThrustersData.xlsx','HiPEP','G28:G37');
    eta3 = xlsread('ThrustersData.xlsx','HiPEP','D28:D37');
        
    
%This is the calculations for when we hold mDot constant    
    pIn1 = powerIn1-min(powerIn1);
    pIn2 = powerIn2-min(powerIn2);
    pIn3 = powerIn3-min(powerIn3);
   
     
    %make a smoother line by using more Pin values
    
    stepSizeB = [(max(pIn1) - min(pIn1))/ INCREMENT_CONSTANT;(max(pIn2) - min(pIn2))/ INCREMENT_CONSTANT;
        (max(pIn3) - min(pIn3))/ INCREMENT_CONSTANT];
    
    %%%sets powerSteps to have the min start at zero for giving it a better
    %%%line fit
    
    powerSteps1 = min(pIn1):stepSizeB(1,1):max(pIn1);
     powerSteps2 = min(pIn2):stepSizeB(2,1):max(pIn2);
      powerSteps3 = min(pIn3):stepSizeB(3,1):max(pIn3);
       
     
    lineFit = polyfit(sqrt(pIn1),eta1,1);
    YMDot1 = polyval(lineFit,sqrt(powerSteps1));
    
    lineFit = polyfit(sqrt(pIn2),eta2,1);
    YMDot2 = polyval(lineFit,sqrt(powerSteps2));
    
    lineFit = polyfit(sqrt(pIn3),eta3,1);
    YMDot3 = polyval(lineFit,sqrt(powerSteps3));
   
    %%%%% set powerSteps back to measued min and max values for ploting
    
    
     powerSteps1 = min(powerIn1):stepSizeB(1,1):max(powerIn1);
      powerSteps2 = min(powerIn2):stepSizeB(2,1):max(powerIn2);
       powerSteps3 = min(powerIn3):stepSizeB(3,1):max(powerIn3);
         
    
    
    
    
%This is the calculations for when we hold uExhaust constant

    %First set with Uex = 90000 m/s
    powerInB1 = xlsread('ThrustersData.xlsx','HiPEP','P3:P4');
    etaB1 = xlsread('ThrustersData.xlsx','HiPEP','M3:M4');
    
    %Second set with Uex = 85000 m/s
    powerInB2 = xlsread('ThrustersData.xlsx','HiPEP','P5:P7');
    etaB2 = xlsread('ThrustersData.xlsx','HiPEP','M5:M7');
   
    %Third set with Uex = 80000 m/s
    powerInB3 = xlsread('ThrustersData.xlsx','HiPEP','P8:P10');
    etaB3 = xlsread('ThrustersData.xlsx','HiPEP','M8:M10');
    
    %Fourth set with Uex = 75000 m/s
    powerInB4 = xlsread('ThrustersData.xlsx','HiPEP','P11:P12');
    etaB4 = xlsread('ThrustersData.xlsx','HiPEP','M11:M12');
    
    
    pInB1 = powerInB1-min(powerInB1);
    pInB2 = powerInB2-min(powerInB2);
    pInB3 = powerInB3-min(powerInB3);
    pInB4 = powerInB4-min(powerInB4);
   
    %make a smoother line by using more Pin values
    
    stepSizeB = [(max(pInB1) - min(pInB1))/ INCREMENT_CONSTANT;(max(pInB2) - min(pInB2))/ INCREMENT_CONSTANT;
        (max(pInB3) - min(pInB4))/ INCREMENT_CONSTANT;(max(pInB4) - min(pInB4))/ INCREMENT_CONSTANT];
    
    %%%sets powerSteps to have the min start at zero for giving it a better
    %%%line fit
    
    powerStepsB1 = min(pInB1):stepSizeB(1,1):max(pInB1);
     powerStepsB2 = min(pInB2):stepSizeB(2,1):max(pInB2);
      powerStepsB3 = min(pInB3):stepSizeB(3,1):max(pInB3);
       powerStepsB4 = min(pInB4):stepSizeB(4,1):max(pInB4);
     
     
    lineFit = polyfit(sqrt(pInB1),etaB1,1);
    YUex1 = polyval(lineFit,sqrt(powerStepsB1));
    
    lineFit = polyfit(sqrt(pInB2),etaB2,1);
    YUex2 = polyval(lineFit,sqrt(powerStepsB2));
    
    lineFit = polyfit(sqrt(pInB3),etaB3,1);
    YUex3 = polyval(lineFit,sqrt(powerStepsB3));
   
    lineFit = polyfit(sqrt(pInB4),etaB4,1);
    YUex4 = polyval(lineFit,sqrt(powerStepsB4));

    %%%%% set powerSteps back to measued min and max values for ploting
    
    
     powerStepsB1 = min(powerInB1):stepSizeB(1,1):max(powerInB1);
      powerStepsB2 = min(powerInB2):stepSizeB(2,1):max(powerInB2);
       powerStepsB3 = min(powerInB3):stepSizeB(3,1):max(powerInB3);
        powerStepsB4 = min(powerInB4):stepSizeB(4,1):max(powerInB4);
         
    
       
    plot(powerSteps1,YMDot1,'c-',powerStepsB1,YUex1,'m-',powerInPlotting,etaPlotting,'k+',powerSteps2,YMDot2,'c-',powerSteps3,YMDot3,'c-',powerStepsB1,YUex1,'m-',powerStepsB2,YUex2,'m-',powerStepsB3,YUex3,'m-',powerStepsB4,YUex4,'m-');
  
    %Labels the mDot Lines
    mDotLabels = {'4 mg/s','6.6 mg/s','7.1 mg/s'};
    text(powerIn1(1,1),YMDot1(1,1),mDotLabels(1,1),'VerticalAlignment','top','HorizontalAlignment','left')
    text(powerIn2(1,1),YMDot2(1,1),mDotLabels(1,2),'VerticalAlignment','top','HorizontalAlignment','left')
    text(powerIn3(1,1),YMDot3(1,1),mDotLabels(1,3),'VerticalAlignment','top','HorizontalAlignment','left')
    
    %Labels the uExhaust Lines
    uExhaustLabels = {'90000 m/s','85000 m/s','80000 m/s','75000 m/s'};
    text(powerInB1(1,1),YUex1(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,1),'VerticalAlignment','top','HorizontalAlignment','right')
    text(powerInB2(1,1),YUex2(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,2),'VerticalAlignment','bottom','HorizontalAlignment','right')
    text(powerInB3(1,1),YUex3(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,3),'VerticalAlignment','bottom','HorizontalAlignment','right')
    text(powerInB4(1,1),YUex4(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,4),'VerticalAlignment','bottom','HorizontalAlignment','right')
    
    title(thrusterName)
    xlabel('Power (W)')
    ylabel('Efficiency')
    axis([9000,38000,0.7,.8])
    legend('Constant mDot','Constant uExhaust','Measured Values','Location','southeast')

    
    
    
    

%%%%%%%%%%%%%%%%%%%% NEXT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif thrusterID == 6
    
    thrusterName = 'NEXT';
    %read in measured point to plot
    powerInPlotting = xlsread('ThrustersData.xlsx','NEXT','A53:A92');
    etaPlotting = xlsread('ThrustersData.xlsx','NEXT','D53:D92');
    
   
    %first set with mDot at 5.76E-06 kg/s
    powerIn1 = xlsread('ThrustersData.xlsx','NEXT', 'G3:G6');
    eta1 = xlsread('ThrustersData.xlsx','NEXT','C3:C6');
    
    
    %Second set with mDot at 5.12E-06 kg/s
     powerIn2 = xlsread('ThrustersData.xlsx','NEXT','G8:G11');
     eta2 = xlsread('ThrustersData.xlsx','NEXT','C8:C11');
   
    
    %Third set with mDot at 4.46E-06 kg/s
    powerIn3 = xlsread('ThrustersData.xlsx','NEXT','G13:G17');
    eta3 = xlsread('ThrustersData.xlsx','NEXT','C13:C17');
   
    
    %Fourth set with mDot at 3.92E-06 kg/s
    powerIn4 = xlsread('ThrustersData.xlsx','NEXT','G19:G23');
    eta4 = xlsread('ThrustersData.xlsx','NEXT','C19:C23');
   
    %Fifth set with mDot at 3.16E-06 kg/s
    powerIn5 = xlsread('ThrustersData.xlsx','NEXT','G25:G29');
    eta5 = xlsread('ThrustersData.xlsx','NEXT','C25:C29');
   
    
    %Sixth set with mDot at 2.60E-06 kg/s
    powerIn6 = xlsread('ThrustersData.xlsx','NEXT','G31:G35');
    eta6 = xlsread('ThrustersData.xlsx','NEXT','C31:C35');
    
    
    %Seventh set with mDot at 2.05E-06 kg/s
    powerIn7 = xlsread('ThrustersData.xlsx','NEXT','G37:G47');
    eta7 = xlsread('ThrustersData.xlsx','NEXT','C37:C47');
 
    
    %Eighth set with mDot at 1.85E-06 kg/s
    powerIn8 = xlsread('ThrustersData.xlsx','NEXT','G49:G49');
    eta8 = xlsread('ThrustersData.xlsx','NEXT','C49:C49');
  
    
    pIn1 = powerIn1-min(powerIn1);
    pIn2 = powerIn2-min(powerIn2);
    pIn3 = powerIn3-min(powerIn3);
    pIn4 = powerIn4-min(powerIn4);
    pIn5 = powerIn5-min(powerIn5);
    pIn6 = powerIn6-min(powerIn6);
    pIn7 = powerIn7-min(powerIn7);
    pIn8 = powerIn8-min(powerIn8);
    
    
     stepSize = [(max(pIn1) - min(pIn1))/ INCREMENT_CONSTANT;(max(pIn2) - min(pIn2))/ INCREMENT_CONSTANT;(max(pIn3) - min(pIn3))/ INCREMENT_CONSTANT;
        (max(pIn4) - min(pIn4))/ INCREMENT_CONSTANT;(max(pIn5) - min(pIn5))/ INCREMENT_CONSTANT;(max(pIn6) - min(pIn6))/ INCREMENT_CONSTANT;
        (max(pIn7) - min(pIn7))/ INCREMENT_CONSTANT;(max(pIn8) - min(pIn8))/ INCREMENT_CONSTANT;];
    
    
    powerSteps1 = max(pIn1):-stepSize(1,1):min(pIn1);
     powerSteps2 = max(pIn2):-stepSize(2,1):min(pIn2);
      powerSteps3 = max(pIn3):-stepSize(3,1):min(pIn3);
       powerSteps4 = max(pIn4):-stepSize(4,1):min(pIn4);
        powerSteps5 = max(pIn5):-stepSize(5,1):min(pIn5);
         powerSteps6 = max(pIn6):-stepSize(6,1):min(pIn6);
          powerSteps7 = max(pIn7):-stepSize(7,1):min(pIn7);
           powerSteps8 = max(pIn8):-stepSize(8,1):min(pIn8);
               
           powerSteps1 = powerSteps1';
            powerSteps2 = powerSteps2';
             powerSteps3 = powerSteps3';
              powerSteps4 = powerSteps4';
               powerSteps5 = powerSteps5';
                powerSteps6 = powerSteps6';
                 powerSteps7 = powerSteps7';
                  powerSteps8 = powerSteps8';
    
     
    
                  
    lineFit = polyfit(sqrt(pIn1),eta1,1);
    YMDot1 = polyval(lineFit,sqrt(powerSteps1));
       
    lineFit = polyfit(sqrt(pIn2),eta2,1);
    YMDot2 = polyval(lineFit,sqrt(powerSteps2));
       
    lineFit = polyfit(sqrt(pIn3),eta3,1);
    YMDot3 = polyval(lineFit,sqrt(powerSteps3));
    
    lineFit = polyfit(sqrt(pIn4),eta4,1);
    YMDot4 = polyval(lineFit,sqrt(powerSteps4));
    
    lineFit = polyfit(sqrt(pIn5),eta5,1);
    YMDot5 = polyval(lineFit,sqrt(powerSteps5));
    
    lineFit = polyfit(sqrt(pIn6),eta6,1);
    YMDot6 = polyval(lineFit,sqrt(powerSteps6));
    
    lineFit = polyfit(sqrt(pIn7),eta7,1);
    YMDot7 = polyval(lineFit,sqrt(powerSteps7));
    
    lineFit = polyfit(sqrt(pIn8),eta8,1);
    YMDot8 = polyval(lineFit,sqrt(powerSteps8));
    
    %%%%%set powerSteps back to measued min and max values for ploting
    
    powerSteps1 = max(powerIn1):-stepSize(1,1):min(powerIn1);
     powerSteps2 = max(powerIn2):-stepSize(2,1):min(powerIn2);
      powerSteps3 = max(powerIn3):-stepSize(3,1):min(powerIn3);
       powerSteps4 = max(powerIn4):-stepSize(4,1):min(powerIn4);
        powerSteps5 = max(powerIn5):-stepSize(5,1):min(powerIn5);
         powerSteps6 = max(powerIn6):-stepSize(6,1):min(powerIn6);
          powerSteps7 = max(powerIn7):-stepSize(7,1):min(powerIn7);
           powerSteps8 = max(powerIn8):-stepSize(8,1):min(powerIn8);
    
   
        
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%'
    %%%%%%%%%%%%%%%% NEXT DATA HOLDING Uex CONSTANT %%%%%%%%%%%%%
    
    
    %First set with Uex = 40000 m/s
    powerInB1 = xlsread('ThrustersData.xlsx','NEXT', 'R3:R8');
    etaB1 = xlsread('ThrustersData.xlsx','NEXT','N3:N8');
    
    %Second  set with Uex = 37500 m/s
    powerInB2 = xlsread('ThrustersData.xlsx','NEXT', 'R9:R15');
    etaB2 = xlsread('ThrustersData.xlsx','NEXT','N9:N15');
    
    %Third  set with Uex = 35000 m/s
    powerInB3 = xlsread('ThrustersData.xlsx','NEXT', 'R16:R22');
    etaB3 = xlsread('ThrustersData.xlsx','NEXT','N16:N22');
    
    %Fourth  set with Uex = 32500 m/s
    powerInB4 = xlsread('ThrustersData.xlsx','NEXT', 'R23:R27');
    etaB4 = xlsread('ThrustersData.xlsx','NEXT','N23:N27');

    
    pInB1 = powerInB1-min(powerInB1);
     pInB2 = powerInB2-min(powerInB2);
      pInB3 = powerInB3-min(powerInB3);
       pInB4 = powerInB4-min(powerInB4);
    
    %make a smoother line by using more Pin values
    
    stepSizeB = [(max(pInB1) - min(pInB1))/ INCREMENT_CONSTANT;(max(pInB2) - min(pInB2))/ INCREMENT_CONSTANT;
        (max(pInB3) - min(pInB3))/ INCREMENT_CONSTANT;(max(pInB4) - min(pInB4))/ INCREMENT_CONSTANT];
    
    %%%sets powerSteps to have the min start at zero for giving it a better
    %%%line fit
    
    powerStepsB1 = min(pInB1):stepSizeB(1,1):max(pInB1);
     powerStepsB2 = min(pInB2):stepSizeB(2,1):max(pInB2);
      powerStepsB3 = min(pInB3):stepSizeB(3,1):max(pInB3);
       powerStepsB4 = min(pInB4):stepSizeB(4,1):max(pInB4);
     
     
    lineFit = polyfit(sqrt(pInB1),etaB1,1);
    YUex1 = polyval(lineFit,sqrt(powerStepsB1));
    
    lineFit = polyfit(sqrt(pInB2),etaB2,1);
    YUex2 = polyval(lineFit,sqrt(powerStepsB2));
    
    lineFit = polyfit(sqrt(pInB3),etaB3,1);
    YUex3 = polyval(lineFit,sqrt(powerStepsB3));
   
    lineFit = polyfit(sqrt(pInB4),etaB4,1);
    YUex4 = polyval(lineFit,sqrt(powerStepsB4));

    %%%%% set powerSteps back to measued min and max values for ploting
    
    
     powerStepsB1 = min(powerInB1):stepSizeB(1,1):max(powerInB1);
      powerStepsB2 = min(powerInB2):stepSizeB(2,1):max(powerInB2);
       powerStepsB3 = min(powerInB3):stepSizeB(3,1):max(powerInB3);
        powerStepsB4 = min(powerInB4):stepSizeB(4,1):max(powerInB4);
         
    
       
    
    plot(powerSteps1,YMDot1,'c-', powerStepsB1, YUex1, '-.r',powerInPlotting,etaPlotting,'k+',powerSteps2,YMDot2,'c-',powerSteps3,YMDot3,'c-',powerSteps4,YMDot4,'c-',powerSteps5,YMDot5,'c-',powerSteps6,YMDot6,'c-',powerSteps7,YMDot7,'c-',powerSteps8,YMDot8,'c-',powerStepsB2, YUex2, '-.r', powerStepsB3, YUex3, '-.r', powerStepsB4, YUex4, '-.r');
    
    mDotLabels = {'5.7 mg/s','5.1 mg/s','4.4 mg/s','3.9 mg/s','3.2 mg/s','2.6 mg/s','2.0 mg/s','1.9 mg/s'};
    text(powerIn1(1,1),YMDot1(1,1),mDotLabels(1,1),'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(powerIn2(1,1),YMDot2(1,1),mDotLabels(1,2),'VerticalAlignment','top','HorizontalAlignment','left')
    text(powerIn3(1,1),YMDot3(1,1),mDotLabels(1,3),'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(powerIn4(1,1),YMDot4(1,1),mDotLabels(1,4),'VerticalAlignment','top','HorizontalAlignment','left')
    text(powerIn5(1,1),YMDot5(1,1),mDotLabels(1,5),'VerticalAlignment','bottom','HorizontalAlignment','right')
    text(powerIn6(1,1),YMDot6(1,1),mDotLabels(1,6),'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(powerIn7(1,1),YMDot7(1,1),mDotLabels(1,7),'VerticalAlignment','bottom','HorizontalAlignment','left')
    %text(powerIn8(1,1),YMDot8(1,1),mDotLabels(1,8),'VerticalAlignment','top','HorizontalAlignment','left')
    
    %Labels the uExhaust Lines
    uExhaustLabels = {'40000 m/s','37500 m/s','35000 m/s','32500 m/s'};
    text(powerInB1(1,1),YUex1(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,1),'VerticalAlignment','top','HorizontalAlignment','left')
    text(powerInB2(1,1),YUex2(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,2),'VerticalAlignment','bottom','HorizontalAlignment','left')
    text(powerInB3(1,1),YUex3(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,3),'VerticalAlignment','top','HorizontalAlignment','left')
    text(powerInB4(1,1),YUex4(1,1+INCREMENT_CONSTANT),uExhaustLabels(1,4),'VerticalAlignment','top','HorizontalAlignment','left')
    
    title(thrusterName)
    xlabel('Power (W)')
    ylabel('Efficiency')
    axis([600,7500,0.3,.8])
    legend('Constant mDot','Constant uExhaust','Measured Values','Location','southeast')
    
    
    
end

end

