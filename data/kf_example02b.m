% kf_example02b.m
%  
%  This Matlab script loads into your work space the problem matrices
%  and data for a linear discrete-time Kalman filter.
%
   Fk     = [  0.81671934103521,  0.08791146849849;...
              -3.47061412053765,  0.70624978972000];     % for all k
   Gammak = [  0.00464254201630;...
               0.08791146849849];                        % for all k
   Hk     = [  2.00000000000000,  0.30000000000000];     % for all k
%
   Qk     =   40.00000000000000;                         % for all k
%  Qk     =    0.40000000000000;      alternate          % for all k
%  Qk     =    0.00400000000000;      alternate          % for all k
   Rk     =    0.00010000000000;                         % for all k
%
   xhat0   = [  0.20000000000000;...
               -2.50000000000000];
   P0      = [  0.25000000000000,  0.08000000000000;...
                0.08000000000000,  0.50000000000000];
%
%  Note that the following arrays define the sample times 
%  [t(1);t(2);t(3);...;t(50)] and the measurements 
%  [z(1);z(2);z(3);...;z(50)]
%
   thist = [1:50]'*0.1;
   zhist = [-0.75352577483108;...
            -1.47835683178354;...
            -1.60356954253130;...
            -1.14319787399192;...
            -0.29770494780324;...
             0.54807996865678;...
             1.07646986214639;...
             1.21029723800666;...
             0.86753373429266;...
             0.26639106635320;...
            -0.39741466924670;...
            -0.81680173977677;...
            -0.89780421626560;...
            -0.60343182449369;...
            -0.17504910931582;...
             0.29126479721498;...
             0.64256054671834;...
             0.69234407038111;...
             0.47952867123181;...
             0.15536903510622;...
            -0.21717489024093;...
            -0.45158274203103;...
            -0.51519554865583;...
            -0.36237832381512;...
            -0.07725284799720;...
             0.19593385632001;...
             0.32312914618026;...
             0.34784096787121;...
             0.23319289376215;...
             0.06293622960977;...
            -0.09870551526524;...
            -0.24060000935609;...
            -0.27989488885854;...
            -0.21677480698919;...
            -0.05506520661488;...
             0.07468611197805;...
             0.16393837131234;...
             0.17032815074154;...
             0.17219253765052;...
             0.10489396937737;...
             0.03254186299003;...
            -0.05701463741298;...
            -0.12708672066528;...
            -0.11839177276129;...
            -0.08759623884611;...
            -0.03288894039672;...
             0.04546176011055;...
             0.10520182441580;...
             0.11771261524308;...
             0.07633914203888];
