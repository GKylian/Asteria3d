#pragma once
#include "MHD.h"
#include <iostream>


bool X_HLLFlux(ld uL[8], ld uR[8], ld SL, ld SR, ld FL[8], ld FR[8], ld *F) {

    //1. Get primitive variables and compute intermediate speed
    //ld wL[8] = { 0 }; if (!toPrimitive(uL, wL)) return false; ld wR[8] = { 0 }; if (!toPrimitive(uR, wR)) return false;
    

    //3. Choose which intermediate flux to return
    if (SL > 0.0) {
        for (int i = 0; i < 8; i++)
            F[i] = FL[i];
        return true;
    }
    else if (SL <= 0.0 && SR >= 0.0) {
        for (int i = 0; i < 8; i++)
            F[i] = (SR*FL[i] - SL*FR[i] + SR*SL*(uR[i]-uL[i])) / (SR-SL);
        return true;
    }
    else if (SR < 0.0) {
        for (int i = 0; i < 8; i++)
            F[i] = FR[i];
        return true;
    }


    return false;
}

bool Y_HLLFlux(ld uL[8], ld uR[8], ld SL, ld SR, ld GL[8], ld GR[8], ld *G) {

    //1. Get primitive variables and compute intermediate speed
    //ld wL[8] = { 0 }; if (!toPrimitive(uL, wL)) return false; ld wR[8] = { 0 }; if (!toPrimitive(uR, wR)) return false;
    

    //3. Choose which intermediate flux to return
    if (SL > 0.0) {
        for (int i = 0; i < 8; i++)
            G[i] = GL[i];
        return true;
    }
    else if (SL <= 0.0 && SR >= 0.0) {
        for (int i = 0; i < 8; i++)
            G[i] = (SR*GL[i] - SL*GR[i] + SR*SL*(uR[i]-uL[i])) / (SR-SL);
        return true;
    }
    else if (SR < 0.0) {
        for (int i = 0; i < 8; i++)
            G[i] = GR[i];
        return true;
    }


    return false;
}


//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------


bool X_HLLCFlux(ld uL[8], ld uR[8], ld SL, ld SR, ld FL[8], ld FR[8], ld *F) {

    //1. Get primitive variables and compute intermediate speed
    ld wL[8] = { 0 }; if(!toPrimitive(uL, wL)) return false; ld wR[8] = { 0 }; if (!toPrimitive(uR, wR)) return false;
    
    ld B2L = uL[5]*uL[5]+uL[6]*uL[6]+uL[7]*uL[7]; ld B2R = uR[5]*uR[5]+uR[6]*uR[6]+uR[7]*uR[7];
    ld PTL = wL[4] + B2L/2.0; ld PTR = wR[4] + B2R/2.0;

    ld SM = ((SR-wR[1])*uR[1] - (SL-wL[1])*uL[1] - PTR + PTL) / ((SR-wR[1])*uR[0] - (SL-wL[1])*uL[0]);
    ld PTsL = PTL + uL[0]*(SL-wL[1])*(SM-wL[1]); ld PTsR = PTR + uR[0]*(SR-wR[1])*(SM-wR[1]);
    ld PTs = (PTsL+PTsR)/2.0;
    //if (PTs != PTR + uR[0]*(SR-wR[1])*(SM-wR[1])) cout << "X_HLLCFlux:: PTs at left and right interfaces not constant." << endl;

    //2. Get intermediate variables
    ld usL[8] = { 0 }; ld usR[8] = { 0 };
    usL[0] = uL[0]*(SL-wL[1])/(SL-SM);      usR[0] = uR[0]*(SR-wR[1])/(SR-SM);
    usL[5] = uL[5]; usR[5] = uR[5];
    usL[1] = SM; usR[1] = SM;

    /*for (int i = 2; i < 8; i++)
    {
        if (i == 5 || i == 4) continue;
        usL[i] = ( SR*uR[i] - SL*uL[i] - FR[i] + FL[i] ) / (SR - SL); usR[i] = usL[i];
    }*/
    if (uL[5] != 0) {


        usL[2] = (  uL[0]*wL[2]*(wL[1]-SL) - uR[0]*wR[2]*(wR[1]-SR) + uL[5]*(uL[6]-uR[6])/2.0  ) / (uL[0]*(wL[1]-SL) - uR[0]*(wR[1]-SR)); usR[2] = usL[2];
        usL[3] = (  uL[0]*wL[3]*(wL[1]-SL) - uR[0]*wR[3]*(wR[1]-SR) + uL[5]*(uL[7]-uR[7])/2.0  ) / (uL[0]*(wL[1]-SL) - uR[0]*(wR[1]-SR)); usR[3] = usL[3];

        usL[6] = (  uL[6]*(wL[1]-SL) - uR[6]*(wR[1]-SR) + uL[5]*(wL[2]-wR[2])  ) / ( SR - SL ); usR[6] = usL[6];
        usL[7] = (  uL[7]*(wL[1]-SL) - uR[7]*(wR[1]-SR) + uL[5]*(wL[3]-wR[3])  ) / ( SR - SL ); usR[7] = usL[7];


        //usL[2] = wL[2] + wL[6]*wL[5]*(SM-wL[1])/(  wL[5]*wL[5] + wL[0]*(wL[1]-SL)*(SL-SM) );
        //usR[2] = wR[2] + wR[6]*wR[5]*(SM-wR[1])/(  wR[5]*wR[5] + wR[0]*(wR[1]-SR)*(SR-SM) );

        //usL[3] = wL[3] + wL[7]*wL[5]*(SM-wL[1])/(  wL[5]*wL[5] + wL[0]*(wL[1]-SL)*(SL-SM) );
        //usR[3] = wR[3] + wR[7]*wR[5]*(SM-wR[1])/(  wR[5]*wR[5] + wR[0]*(wR[1]-SR)*(SR-SM) );

        //usL[6] = wL[6]*( wL[5]*wL[5] - wL[0]*powl(wL[1]-SL,2) ) / (  wL[5]*wL[5] + wL[0]*(wL[1]-SL)*(SL-SM) );
        //usR[6] = wR[6]*( wR[5]*wR[5] - wR[0]*powl(wR[1]-SR,2) ) / (  wR[5]*wR[5] + wR[0]*(wR[1]-SR)*(SR-SM) );

        //usL[7] = wL[7]*( wL[5]*wL[5] - wL[0]*powl(wL[1]-SL,2) ) / (  wL[5]*wL[5] + wL[0]*(wL[1]-SL)*(SL-SM) );
        //usR[7] = wR[7]*( wR[5]*wR[5] - wR[0]*powl(wR[1]-SR,2) ) / (  wR[5]*wR[5] + wR[0]*(wR[1]-SR)*(SR-SM) );


        //usL[2] = wL[2] + uL[5]*(uL[6]-usL[6])/(uL[0]*(SL-wL[1])); usR[2] = wR[2] + uR[5]*(uR[6]-usR[6])/(uR[0]*(SR-wR[1]));//   This correction doesn't change anything in RJ2a.
        //usL[3] = wL[3] + uL[5]*(uL[7]-usL[7])/(uL[0]*(SL-wL[1])); usR[3] = wR[3] + uR[5]*(uR[7]-usR[7])/(uR[0]*(SR-wR[1]));


        usL[1] *= usL[0]; usR[1] *= usR[0]; //Velocity to momentum
        usL[2] *= usL[0]; usR[2] *= usR[0];
        usL[3] *= usL[0]; usR[3] *= usR[0];

        ld uL_dot_BL = wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]; ld uR_dot_BR = wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7];
        ld usL_dot_BsL = usL[1]/usL[0]*usL[5] + usL[2]/usL[0]*usL[6] +usL[3]/usL[0]*usL[7];
        ld usR_dot_BsR = usR[1]/usR[0]*usR[5] + usR[2]/usR[0]*usR[6] +usR[3]/usR[0]*usR[7];

        usL[4] = ((SL-wL[1])*uL[4] - PTL*wL[1] + PTs*SM + usL[5]*(uL_dot_BL-usL_dot_BsL)) / (SL-SM);   usR[4] = ((SR-wR[1])*uR[4] - PTR*wR[1] + PTs*SM + usR[5]*(uR_dot_BR-usR_dot_BsR)) / (SR-SM);


    }
    else {

        usL[2] = wL[2]; usR[2] = wR[2];
        usL[3] = wL[3]; usR[3] = wR[3];

        usL[6] = uL[6]*(SL-wL[1])/(SL-SM);      usR[6] = uR[6]*(SR-wR[1])/(SR-SM);
        usL[7] = uL[7]*(SL-wL[1])/(SL-SM);      usR[7] = uR[7]*(SR-wR[1])/(SR-SM);

        //usL[2] = wL[2] + uL[5]*(uL[6]-usL[6])/(uL[0]*(SL-wL[1])); usR[2] = wR[2] + uR[5]*(uR[6]-usR[6])/(uR[0]*(SR-wR[1]));   //The correction fails.
        //usL[3] = wL[3] + uL[5]*(uL[7]-usL[7])/(uL[0]*(SL-wL[1])); usR[3] = wR[3] + uR[5]*(uR[7]-usR[7])/(uR[0]*(SR-wR[1]));


        usL[1] *= usL[0]; usR[1] *= usR[0]; //Velocity to momentum
        usL[2] *= usL[0]; usR[2] *= usR[0];
        usL[3] *= usL[0]; usR[3] *= usR[0];


        /*usL[4] = ((SL-wL[1])*uL[4] - PTL*wL[1]) / (SL-SM);   usR[4] = ((SR-wR[1])*uR[4] - PTR*wR[1]) / (SR-SM);*/

        ld uL_dot_BL = wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]; ld uR_dot_BR = wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7];
        ld usL_dot_BsL = usL[1]/usL[0]*usL[5] + usL[2]/usL[0]*usL[6] +usL[3]/usL[0]*usL[7];
        ld usR_dot_BsR = usR[1]/usR[0]*usR[5] + usR[2]/usR[0]*usR[6] +usR[3]/usR[0]*usR[7];

        usL[4] = ((SL-wL[1])*uL[4] - PTL*wL[1] + PTs*SM + usL[5]*(uL_dot_BL-usL_dot_BsL)) / (SL-SM);   usR[4] = ((SR-wR[1])*uR[4] - PTR*wR[1] + PTs*SM + usR[5]*(uR_dot_BR-usR_dot_BsR)) / (SR-SM);


    }

    
    

    //3. Choose which intermediate flux to return
    if (SL > 0) {
        for (int i = 0; i < 8; i++)
            F[i] = FL[i];
    }
    else if (SL <= 0 && SM >= 0) {
        for (int i = 0; i < 8; i++)
            F[i] = (SR*FL[i] - SL*FR[i] + SR*SL*(uR[i]-uL[i])) / (SR-SL)  -  SL*(SR-SM)/(SR-SL)*(usR[i]-usL[i]);
            //F[i] = FL[i] + SL*(usL[i] - uL[i]);
            //F[i] = FsL[i];
    }
    else if (SM <= 0 && SR >= 0) {
        for (int i = 0; i < 8; i++)
            F[i] = (SR*FL[i] - SL*FR[i] + SR*SL*(uR[i]-uL[i])) / (SR-SL)  +  SR*(SM-SL)/(SR-SL)*(usR[i]-usL[i]);
            //F[i] = FR[i] + SR*(usR[i] - uR[i]);
            //F[i] = FsR[i];
    }
    else if (SR < 0) {
        for (int i = 0; i < 8; i++)
            F[i] = FR[i];
    }


    return true;
}



bool Y_HLLCFlux(ld uL[8], ld uR[8], ld SL, ld SR, ld GL[8], ld GR[8], ld *G) {

    //1. Get primitive variables and compute intermediate speed
    ld wL[8] = { 0 }; if(!toPrimitive(uL, wL)) return false; ld wR[8] = { 0 }; if (!toPrimitive(uR, wR)) return false;
    
    ld B2L = uL[5]*uL[5]+uL[6]*uL[6]+uL[7]*uL[7]; ld B2R = uR[5]*uR[5]+uR[6]*uR[6]+uR[7]*uR[7];
    ld PTL = wL[4] + B2L/2.0; ld PTR = wR[4] + B2R/2.0;

    ld SM = ((SR-wR[2])*uR[2] - (SL-wL[2])*uL[2] - PTR + PTL) / ((SR-wR[2])*uR[0] - (SL-wL[2])*uL[0]);
    ld PTsL = PTL + uL[0]*(SL-wL[2])*(SM-wL[2]); ld PTsR = PTR + uR[0]*(SR-wR[2])*(SM-wR[2]);
    ld PTs = (PTsL+PTsR)/2.0;
    //if (PTs != PTR + uR[0]*(SR-wR[1])*(SM-wR[1])) cout << "X_HLLCGlux:: PTs at left and right interfaces not constant." << endl;

    //2. Get intermediate variables
    ld usL[8] = { 0 }; ld usR[8] = { 0 };
    usL[0] = uL[0]*(SL-wL[2])/(SL-SM);      usR[0] = uR[0]*(SR-wR[2])/(SR-SM);
    usL[6] = uL[6]; usR[6] = uR[6];
    usL[2] = SM; usR[2] = SM;

    /*for (int i = 2; i < 8; i++)
    {
        if (i == 5 || i == 4) continue;
        usL[i] = ( SR*uR[i] - SL*uL[i] - GR[i] + GL[i] ) / (SR - SL); usR[i] = usL[i];
    }*/
    if (uL[6] != 0) {
        usL[1] = (  uL[0]*wL[1]*(wL[2]-SL) - uR[0]*wR[1]*(wR[2]-SR) + uL[6]*(uL[5]-uR[5])/2.0  ) / (uL[0]*(wL[2]-SL) - uR[0]*(wR[2]-SR)); usR[1] = usL[1];
        usL[3] = (  uL[0]*wL[3]*(wL[2]-SL) - uR[0]*wR[3]*(wR[2]-SR) + uL[6]*(uL[7]-uR[7])/2.0  ) / (uL[0]*(wL[2]-SL) - uR[0]*(wR[2]-SR)); usR[3] = usL[3];

        usL[5] = (  uL[5]*(wL[2]-SL) - uR[5]*(wR[2]-SR) + uL[6]*(wL[1]-wR[1])  ) / ( SR - SL ); usR[5] = usL[5];
        usL[7] = (  uL[7]*(wL[2]-SL) - uR[7]*(wR[2]-SR) + uL[6]*(wL[3]-wR[3])  ) / ( SR - SL ); usR[7] = usL[7];


        /*usL[1] = wL[1] + wL[5]*wL[6]*(SM-wL[2])/(wL[6]*wL[6] + wL[0]*(wL[2]-SL)*(SL-SM));
        usR[1] = wR[1] + wR[5]*wR[6]*(SM-wR[2])/(wR[6]*wR[6] + wR[0]*(wR[2]-SR)*(SR-SM));
        
        usL[3] = wL[3] + wL[7]*wL[6]*(SM-wL[2])/(wL[6]*wL[6] + wL[0]*(wL[2]-SL)*(SL-SM));
        usR[3] = wR[3] + wR[7]*wR[6]*(SM-wR[2])/(wR[6]*wR[6] + wR[0]*(wR[2]-SR)*(SR-SM));

        usL[5] = wL[5]*(wL[6]*wL[6] - wL[0]*powl(wL[2]-SL, 2)) / (wL[6]*wL[6] + wL[0]*(wL[2]-SL)*(SL-SM));
        usR[5] = wR[5]*(wR[6]*wR[6] - wR[0]*powl(wR[2]-SR, 2)) / (wR[6]*wR[6] + wR[0]*(wR[2]-SR)*(SR-SM));

        usL[7] = wL[7]*(wL[6]*wL[6] - wL[0]*powl(wL[2]-SL, 2)) / (wL[6]*wL[6] + wL[0]*(wL[2]-SL)*(SL-SM));
        usR[7] = wR[7]*(wR[6]*wR[6] - wR[0]*powl(wR[2]-SR, 2)) / (wR[6]*wR[6] + wR[0]*(wR[2]-SR)*(SR-SM));*/


        //usL[1] = wL[1] + uL[6]*(uL[5]-usL[5])/(uL[0]*(SL-wL[2])); usR[1] = wR[1] + uR[6]*(uR[5]-usR[5])/(uR[0]*(SR-wR[2]));//   This correction doesn't change anything in RJ2a.
        //usL[3] = wL[3] + uL[6]*(uL[7]-usL[7])/(uL[0]*(SL-wL[2])); usR[3] = wR[3] + uR[6]*(uR[7]-usR[7])/(uR[0]*(SR-wR[2]));



        usL[1] *= usL[0]; usR[1] *= usR[0]; //Velocity to momentum
        usL[2] *= usL[0]; usR[2] *= usR[0];
        usL[3] *= usL[0]; usR[3] *= usR[0];

        ld uL_dot_BL = wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]; ld uR_dot_BR = wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7];
        ld usL_dot_BsL = usL[1]/usL[0]*usL[5] + usL[2]/usL[0]*usL[6] +usL[3]/usL[0]*usL[7];
        ld usR_dot_BsR = usR[1]/usR[0]*usR[5] + usR[2]/usR[0]*usR[6] +usR[3]/usR[0]*usR[7];

        usL[4] = ((SL-wL[2])*uL[4] - PTL*wL[2] + PTs*SM + usL[6]*(uL_dot_BL-usL_dot_BsL)) / (SL-SM);   usR[4] = ((SR-wR[2])*uR[4] - PTR*wR[2] + PTs*SM + usR[6]*(uR_dot_BR-usR_dot_BsR)) / (SR-SM);


    }
    else {

        usL[1] = wL[1]; usR[1] = wR[1];
        usL[3] = wL[3]; usR[3] = wR[3];

        usL[5] = uL[5]*(SL-wL[2])/(SL-SM);      usR[5] = uR[5]*(SR-wR[2])/(SR-SM);
        usL[7] = uL[7]*(SL-wL[2])/(SL-SM);      usR[7] = uR[7]*(SR-wR[2])/(SR-SM);


        




        usL[1] *= usL[0]; usR[1] *= usR[0]; //Velocity to momentum
        usL[2] *= usL[0]; usR[2] *= usR[0];
        usL[3] *= usL[0]; usR[3] *= usR[0];


        /*usL[4] = ((SL-wL[1])*uL[4] - PTL*wL[1]) / (SL-SM);   usR[4] = ((SR-wR[1])*uR[4] - PTR*wR[1]) / (SR-SM);*/

        ld uL_dot_BL = wL[1]*wL[5]+wL[2]*wL[6]+wL[3]*wL[7]; ld uR_dot_BR = wR[1]*wR[5]+wR[2]*wR[6]+wR[3]*wR[7];
        ld usL_dot_BsL = usL[1]/usL[0]*usL[5] + usL[2]/usL[0]*usL[6] +usL[3]/usL[0]*usL[7];
        ld usR_dot_BsR = usR[1]/usR[0]*usR[5] + usR[2]/usR[0]*usR[6] +usR[3]/usR[0]*usR[7];

        usL[4] = ((SL-wL[2])*uL[4] - PTL*wL[2] + PTs*SM + usL[6]*(uL_dot_BL-usL_dot_BsL)) / (SL-SM);   usR[4] = ((SR-wR[2])*uR[4] - PTR*wR[2] + PTs*SM + usR[6]*(uR_dot_BR-usR_dot_BsR)) / (SR-SM);


    }

    
    

    //3. Choose which intermediate flux to return
    if (SL > 0) {
        for (int i = 0; i < 8; i++)
            G[i] = GL[i];
    }
    else if (SL <= 0 && SM >= 0) {
        for (int i = 0; i < 8; i++)
            G[i] = (SR*GL[i] - SL*GR[i] + SR*SL*(uR[i]-uL[i])) / (SR-SL)  -  SL*(SR-SM)/(SR-SL)*(usR[i]-usL[i]);
            //G[i] = GL[i] + SL*(usL[i] - uL[i]);
            //G[i] = GsL[i];
    }
    else if (SM <= 0 && SR >= 0) {
        for (int i = 0; i < 8; i++)
            G[i] = (SR*GL[i] - SL*GR[i] + SR*SL*(uR[i]-uL[i])) / (SR-SL)  +  SR*(SM-SL)/(SR-SL)*(usR[i]-usL[i]);
            //G[i] = GR[i] + SR*(usR[i] - uR[i]);
            //G[i] = GsR[i];
    }
    else if (SR < 0) {
        for (int i = 0; i < 8; i++)
            G[i] = GR[i];
    }


    return true;
}





//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------
//TODO:---------------------------------------------------------------------------------------------------------------------------------------------------------






bool X_HLLDFlux(ld uL[8], ld uR[8], ld SL, ld SR, ld FL[8], ld FR[8], ld *F) {

    if (uL[5] != uR[5]) {
        cout << "X_HLLDFlux:: Bx is different on the left and on the right !!!" << endl; return false;
    }
    ld Bx = uL[5];

    //1. Get primitive variables and compute intermediate speed
    ld wL[8] = { 0 }; if (!toPrimitive(uL, wL)) return false; ld wR[8] = { 0 }; if (!toPrimitive(uR, wR)) return false;
    ld PTL = wL[4] + (uL[5]*uL[5]+uL[6]*uL[6]+uL[7]*uL[7])/2.0; ld PTR = wR[4] + (uR[5]*uR[5]+uR[6]*uR[6]+uR[7]*uR[7])/2.0;

    ld SM = ((SR-wR[1])*uR[1] - (SL-wL[1])*uL[1] - PTR + PTL) / ((SR-wR[1])*uR[0] - (SL-wL[1])*uL[0]);

    //2. Get the exterior intermediate flux (usL=u*L and usR=u*R)
    ld usL[8] = { 0 }; ld usR[8] = { 0 };
    ld PsT = (  (SR-wR[1])*uR[0]*PTL - (SL-wL[1])*uL[0]*PTR + uL[0]*uR[0]*(SR-wR[1])*(SL-wL[1])*(wR[1]-wL[1])  ) / (  (SR-wR[1])*uR[0] - (SL-wL[1])*uL[0]  );


    usL[0] = uL[0]*(SL-wL[1])/(SL-SM);   usR[0] = uR[0]*(SR-wR[1])/(SR-SM);

    usL[1] = SM*usL[0];   usR[1] = SM*usR[0];

    usL[2] = wL[2] - Bx*uL[6] * (SM-wL[1]) / ( uL[0]*(SL-wL[1])*(SL-SM) - Bx*Bx );   usR[2] = wR[2] - Bx*uR[6] * (SM-wR[1]) / ( uR[0]*(SR-wR[1])*(SR-SM) - Bx*Bx );
    usL[2] *= usL[0];   usR[2] *= usR[0]; //Velocity to momentum

    usL[3] = wL[3] - Bx*uL[7] * (SM-wL[1]) / ( uL[0]*(SL-wL[1])*(SL-SM) - Bx*Bx );   usR[3] = wR[3] - Bx*uR[7] * (SM-wR[1]) / ( uR[0]*(SR-wR[1])*(SR-SM) - Bx*Bx );
    usL[3] *= usL[0];   usR[3] *= usR[0]; //Velocity to momentum

    usL[5] = Bx;   usR[5] = Bx;

    usL[6] = uL[6] * ( uL[0]*powl(SL-wL[1], 2) - Bx*Bx ) / ( uL[0]*(SL-wL[1])*(SL-SM) - Bx*Bx );   usR[6] = uR[6] * ( uR[0]*powl(SR-wR[1], 2) - Bx*Bx ) / ( uR[0]*(SR-wR[1])*(SR-SM) - Bx*Bx );

    usL[7] = uL[7] * ( uL[0]*powl(SL-wL[1], 2) - Bx*Bx ) / ( uL[0]*(SL-wL[1])*(SL-SM) - Bx*Bx );   usR[7] = uR[7] * ( uR[0]*powl(SR-wR[1], 2) - Bx*Bx ) / ( uR[0]*(SR-wR[1])*(SR-SM) - Bx*Bx );

    /*cout << "ByL: " << uL[6] << " --> " << usL[6] << endl;
    cout << "ByR: " << uR[6] << " --> " << usR[6] << endl;
    cout << "BzL: " << uL[7] << " --> " << usL[7] << endl;
    cout << "BzR: " << uR[7] << " --> " << usR[7] << endl;*/

    //Avoid NaN:
    if (isnan(usL[2]) || isnan(usL[3]) || isnan(usL[6]) || isnan(usL[7])) {
        //TODO: is it *usL[0] or *uL[0] ?
        //cout << "X-HLLD:One is NaN" << endl;
        usL[2] = wL[2]*usL[0];   usR[2] = wR[2]*usR[0];
        usL[3] = wL[3]*usL[0];   usR[3] = wR[3]*usR[0];
        /*usL[6] = uL[6];   usR[6] = uR[6];
        usL[7] = uL[7];   usR[7] = uR[7];*/
        usL[6] = 0.0;   usR[6] = 0.0;
        usL[7] = 0.0;   usR[7] = 0.0;

    }

    ld uL_dot_BL = wL[1]*Bx+wL[2]*wL[6]+wL[3]*wL[7]; ld uR_dot_BR = wR[1]*Bx+wR[2]*wR[6]+wR[3]*wR[7];
    ld usL_dot_BsL = usL[1]/usL[0]*Bx + usL[2]/usL[0]*usL[6] +usL[3]/usL[0]*usL[7];
    ld usR_dot_BsR = usR[1]/usR[0]*Bx + usR[2]/usR[0]*usR[6] +usR[3]/usR[0]*usR[7];

    usL[4] = (  (SL-wL[1])*uL[4] - PTL*wL[1] + PsT*SM + Bx*(uL_dot_BL-usL_dot_BsL)  ) / (SL-SM);   usR[4] = (  (SR-wR[1])*uR[4] - PTR*wR[1] + PsT*SM + Bx*(uR_dot_BR-usR_dot_BsR)  ) / (SR-SM);


    //If Bx == 0, there's no need to go further (SsL=SsR=SM)
    if (Bx == 0) {
        if (SL > 0) {
            for (int i = 0; i < 8; i++)
                F[i] = FL[i];
            for (int i = 0; i < 8; i++) {
                if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
                //if (fabsl(F[i]) < 1e-12) F[i] = 0;
                
            }
            return true;
        }
        if (SL <= 0 && SM >= 0) {
            ld FsL[8] = { 0 };
            for (int i = 0; i < 8; i++)
                F[i] = FL[i] + SL*(usL[i]-uL[i]);
            for (int i = 0; i < 8; i++) {
                if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
                //if (fabsl(F[i]) < 1e-12) F[i] = 0;
                
            }
            return true;
        }

        if (SR >= 0 && SM <= 0) {
            for (int i = 0; i < 8; i++)
                F[i] = FR[i] + SR*(usR[i]-uR[i]);
            for (int i = 0; i < 8; i++) {
                if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
                //if (fabsl(F[i]) < 1e-12) F[i] = 0;
                
            }
            return true;
        }
        if (SR < 0) {
            for (int i = 0; i < 8; i++)
                F[i] = FR[i];
            for (int i = 0; i < 8; i++) {
                if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
                //if (fabsl(F[i]) < 1e-12) F[i] = 0;
                
            }
            return true;
        }
    }



    //2. Compute the intermediate wavespeeds SsL and SsR:
    ld SsL = SM - fabsl(Bx)/sqrtl(usL[0]);   ld SsR = SM + fabsl(Bx)/sqrtl(usR[0]);


    //3. Compute the second intermediate fluxes (ussL=u**L and ussR=u**R)
    ld ussL[8] = { 0 }; ld ussR[8] = { 0 };
    ld rtrhoL = sqrtl(usL[0]); ld rtrhoR = sqrtl(usR[0]);

    ussL[0] = usL[0]; ussR[0] = usR[0];

    ussL[1] = usL[1]; ussR[1] = usR[1];

    ussL[2] = ( rtrhoL*usL[2]/usL[0] + rtrhoR*usR[2]/usR[0] + (usR[6]-usL[6])*sgn(Bx) ) / (rtrhoL + rtrhoR); ussR[2] = ussL[2];
    ussL[2] *= ussL[0]; ussR[2] *= ussR[0];

    ussL[3] = ( rtrhoL*usL[3]/usL[0] + rtrhoR*usR[3]/usR[0] + (usR[7]-usL[7])*sgn(Bx) ) / (rtrhoL + rtrhoR); ussR[3] = ussL[3];
    ussL[3] *= ussL[0]; ussR[3] *= ussR[0];

    ussL[5] = Bx; ussR[5] = Bx;

    ussL[6] = ( rtrhoL*usR[6] + rtrhoR*usL[6] + rtrhoL*rtrhoR*(usR[2]/usR[0]-usL[2]/usL[0])*sgn(Bx) ) / (rtrhoL + rtrhoR); ussR[6] = ussL[6];

    ussL[7] = ( rtrhoL*usR[7] + rtrhoR*usL[7] + rtrhoL*rtrhoR*(usR[3]/usR[0]-usL[3]/usL[0])*sgn(Bx) ) / (rtrhoL + rtrhoR); ussR[7] = ussL[7];

    usL_dot_BsL = usL[1]/usL[0]*Bx + usL[2]/usL[0]*usL[6] +usL[3]/usL[0]*usL[7];
    usR_dot_BsR = usR[1]/usR[0]*Bx + usR[2]/usR[0]*usR[6] +usR[3]/usR[0]*usR[7];
    ld uss_dot_Bss = ussL[1]/ussL[0]*ussL[5] + ussL[2]/ussL[0]*ussL[6] + ussL[3]/ussL[0]*ussL[7];

    ussL[4] = usL[4] - rtrhoL*(usL_dot_BsL-uss_dot_Bss)*sgn(Bx);   ussL[4] = usR[4] + rtrhoR*(usR_dot_BsR-uss_dot_Bss)*sgn(Bx);

    if (SL > 0) {
        for (int i = 0; i < 8; i++)
            F[i] = FL[i];
        for (int i = 0; i < 8; i++) {
            if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
            //if (fabsl(F[i]) < 1e-12) F[i] = 0;
        }
        return true;
    }
    if (SL <= 0 && SsL >= 0) {
        ld FsL[8] = { 0 };
        for (int i = 0; i < 8; i++)
            F[i] = FL[i] + SL*(usL[i]-uL[i]);
        for (int i = 0; i < 8; i++) {
            if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
            //if (fabsl(F[i]) < 1e-12) F[i] = 0;
            
        }
        return true;
    }
    if (SM >= 0 && SsL <= 0) {
        for (int i = 0; i < 8; i++)
            F[i] = FL[i] + SL*(usL[i]-uL[i])  +  SsL*(ussL[i]-usL[i]);
        for (int i = 0; i < 8; i++) {
            if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
            //if (fabsl(F[i]) < 1e-12) F[i] = 0;
            
        }
        return true;
    }

    if (SM <= 0 && SsR >= 0) {
        for (int i = 0; i < 8; i++)
            F[i] = FR[i] + SR*(usR[i]-uR[i])  +  SsR*(ussR[i]-usR[i]);
        for (int i = 0; i < 8; i++) {
            if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
            //if (fabsl(F[i]) < 1e-12) F[i] = 0;
            
        }
        return true;
    }
    if (SR >= 0 && SsR <= 0) {
        for (int i = 0; i < 8; i++)
            F[i] = FR[i] + SR*(usR[i]-uR[i]);
        for (int i = 0; i < 8; i++) {
            if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
            //if (fabsl(F[i]) < 1e-12) F[i] = 0;
            
        }
        return true;
    }
    if (SR < 0) {
        for (int i = 0; i < 8; i++)
            F[i] = FR[i];
        for (int i = 0; i < 8; i++) {
            if (isnan(F[i])) { cout << "X-HLLD: F[" << i << "] is NaN !"; return false; }
            //if (fabsl(F[i]) < 1e-12) F[i] = 0;
            
        }
        return true;
    }


    return true;
}


bool Y_HLLDFlux(ld uL[8], ld uR[8], ld SL, ld SR, ld GL[8], ld GR[8], ld *G) {

    if (uL[6] != uR[6]) {
        cout << "Y_HLLDFlux:: By is different on the left and on the right !!!" << endl; return false;
    }
    ld By = uL[6];

    //1. Get primitive variables and compute intermediate speed
    ld wL[8] = { 0 }; if (!toPrimitive(uL, wL)) return false; ld wR[8] = { 0 }; if (!toPrimitive(uR, wR)) return false;
    ld PTL = wL[4] + (uL[5]*uL[5]+uL[6]*uL[6]+uL[7]*uL[7])/2.0; ld PTR = wR[4] + (uR[5]*uR[5]+uR[6]*uR[6]+uR[7]*uR[7])/2.0;

    ld SM = ((SR-wR[2])*uR[2] - (SL-wL[2])*uL[2] - PTR + PTL) / ((SR-wR[2])*uR[0] - (SL-wL[2])*uL[0]);

    //2. Get the exterior intermediate flux (usL=u*L and usR=u*R)
    ld usL[8] = { 0 }; ld usR[8] = { 0 };
    ld PsT = ((SR-wR[2])*uR[0]*PTL - (SL-wL[2])*uL[0]*PTR + uL[0]*uR[0]*(SR-wR[2])*(SL-wL[2])*(wR[2]-wL[2])) / ((SR-wR[2])*uR[0] - (SL-wL[2])*uL[0]);


    usL[0] = uL[0]*(SL-wL[2])/(SL-SM);   usR[0] = uR[0]*(SR-wR[2])/(SR-SM);

    usL[1] = wL[1] - By*uL[5] * (SM-wL[2]) / (uL[0]*(SL-wL[2])*(SL-SM) - By*By);   usR[1] = wR[1] - By*uR[5] * (SM-wR[2]) / (uR[0]*(SR-wR[2])*(SR-SM) - By*By);
    usL[1] *= usL[0];   usR[1] *= usR[0]; //Velocity to momentum

    usL[2] = SM*usL[0];   usR[2] = SM*usR[0];

    usL[3] = wL[3] - By*uL[7] * (SM-wL[2]) / (uL[0]*(SL-wL[2])*(SL-SM) - By*By);   usR[3] = wR[3] - By*uR[7] * (SM-wR[2]) / (uR[0]*(SR-wR[2])*(SR-SM) - By*By);
    usL[3] *= usL[0];   usR[3] *= usR[0]; //Velocity to momentum

    usL[5] = uL[5] * (uL[0]*powl(SL-wL[2], 2) - By*By) / (uL[0]*(SL-wL[2])*(SL-SM) - By*By);   usR[5] = uR[5] * (uR[0]*powl(SR-wR[2], 2) - By*By) / (uR[0]*(SR-wR[2])*(SR-SM) - By*By);

    usL[6] = By;   usR[6] = By;

    usL[7] = uL[7] * (uL[0]*powl(SL-wL[2], 2) - By*By) / (uL[0]*(SL-wL[2])*(SL-SM) - By*By);   usR[7] = uR[7] * (uR[0]*powl(SR-wR[2], 2) - By*By) / (uR[0]*(SR-wR[2])*(SR-SM) - By*By);

    //Avoid NaN:
    if (isnan(usL[2]) || isnan(usL[3]) || isnan(usL[6]) || isnan(usL[7])) {
        //TODO: is it *usL[0] or *uL[0] ?
        //cout << "Y-HLLD:One is NaN" << endl;
        usL[1] = wL[1]*usL[0];   usR[1] = wR[1]*usR[0];
        usL[3] = wL[3]*usL[0];   usR[3] = wR[3]*usR[0];
        usL[5] = uL[5];   usR[5] = uR[5];
        usL[7] = uL[7];   usR[7] = uR[7];
    }

    ld uL_dot_BL = wL[1]*wL[5]+wL[2]*By+wL[3]*wL[7]; ld uR_dot_BR = wR[1]*wR[5]+wR[2]*By+wR[3]*wR[7];
    ld usL_dot_BsL = usL[1]/usL[0]*usL[5] + usL[2]/usL[0]*By +usL[3]/usL[0]*usL[7];
    ld usR_dot_BsR = usR[1]/usR[0]*usR[5] + usR[2]/usR[0]*By +usR[3]/usR[0]*usR[7];

    usL[4] = ((SL-wL[2])*uL[4] - PTL*wL[2] + PsT*SM + By*(uL_dot_BL-usL_dot_BsL)) / (SL-SM);   usR[4] = ((SR-wR[2])*uR[4] - PTR*wR[2] + PsT*SM + By*(uR_dot_BR-usR_dot_BsR)) / (SR-SM);


    //If By=0, there's no need to go further (SsL=SsR=SM)
    if (By==0) {
        if (SL > 0) {
            for (int i = 0; i < 8; i++)
                G[i] = GL[i];
            for (int i = 0; i < 8; i++) {
                if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
                //if (fabsl(G[i]) < 1e-12) G[i] = 0;

                /*if (GL[i] == 0 && GR[i] == 0 && G[i] != 0)
                    cout << "1. G[" << i << "] should be zero, but is " << G[i] << " but GL = " << GL[i] << " and GR = " << GR[i]  << endl;*/
            }
            return true;
        }
        if (SL <= 0 && SM >= 0) {
            for (int i = 0; i < 8; i++)
                G[i] = GL[i] + SL*(usL[i]-uL[i]);
            for (int i = 0; i < 8; i++) {
                if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
                //if (fabsl(G[i]) < 1e-12) G[i] = 0;

                /*if (GL[i] == 0 && GR[i] == 0 && G[i] != 0)
                    cout << "2. G[" << i << "] should be zero, but is " << G[i] << " but GL = " << GL[i] << " and GR = " << GR[i] << endl;*/
            }
            return true;
        }

        if (SR >= 0 && SM <= 0) {
            for (int i = 0; i < 8; i++)
                G[i] = GR[i] + SR*(usR[i]-uR[i]);
            for (int i = 0; i < 8; i++) {
                if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
                //if (fabsl(G[i]) < 1e-12) G[i] = 0;

                /*if (GL[i] == 0 && GR[i] == 0 && G[i] != 0)
                    cout << "3. G[" << i << "] should be zero, but is " << G[i] << " but GL = " << usR[i] << " and GR = " << uR[i]  << endl;*/
            }
            return true;
        }
        if (SR < 0) {
            for (int i = 0; i < 8; i++)
                G[i] = GR[i];
            for (int i = 0; i < 8; i++) {
                if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
                //if (fabsl(G[i]) < 1e-12) G[i] = 0;

                /*if (GL[i] == 0 && GR[i] == 0 && G[i] != 0)
                    cout << "4. G[" << i << "] should be zero, but is " << G[i] << " but GL = " << GL[i] << " and GR = " << GR[i]  << endl;*/
            }
            return true;
        }
    }



    //2. Compute the intermediate wavespeeds SsL and SsR:
    ld SsL = SM - fabsl(By)/sqrtl(usL[0]);   ld SsR = SM + fabsl(By)/sqrtl(usR[0]);


    //3. Compute the second intermediate fluxes (ussL=u**L and ussR=u**R)
    ld ussL[8] = { 0 }; ld ussR[8] = { 0 };
    ld rtrhoL = sqrtl(usL[0]); ld rtrhoR = sqrtl(usR[0]);

    ussL[0] = usL[0]; ussR[0] = usR[0];

    ussL[1] = (rtrhoL*usL[1]/usL[0] + rtrhoR*usR[1]/usR[0] + (usR[5]-usL[5])*sgn(By)) / (rtrhoL + rtrhoR); ussR[1] = ussL[1];
    ussL[1] *= ussL[0]; ussR[1] *= ussR[0];

    ussL[2] = usL[2]; ussR[2] = usR[2];

    ussL[3] = (rtrhoL*usL[3]/usL[0] + rtrhoR*usR[3]/usR[0] + (usR[7]-usL[7])*sgn(By)) / (rtrhoL + rtrhoR); ussR[3] = ussL[3];
    ussL[3] *= ussL[0]; ussR[3] *= ussR[0];

    ussL[5] = (rtrhoL*usR[5] + rtrhoR*usL[5] + rtrhoL*rtrhoR*(usR[1]/usR[0]-usL[1]/usL[0])*sgn(By)) / (rtrhoL + rtrhoR); ussR[5] = ussL[5];

    ussL[6] = By; ussR[6] = By;

    ussL[7] = (rtrhoL*usR[7] + rtrhoR*usL[7] + rtrhoL*rtrhoR*(usR[3]/usR[0]-usL[3]/usL[0])*sgn(By)) / (rtrhoL + rtrhoR); ussR[7] = ussL[7];

    usL_dot_BsL = usL[1]/usL[0]*usL[5] + usL[2]/usL[0]*By +usL[3]/usL[0]*usL[7];
    usR_dot_BsR = usR[1]/usR[0]*usR[5] + usR[2]/usR[0]*By +usR[3]/usR[0]*usR[7];
    ld uss_dot_Bss = ussL[1]/ussL[0]*ussL[5] + ussL[2]/ussL[0]*ussL[6] + ussL[3]/ussL[0]*ussL[7];

    ussL[4] = usL[4] - rtrhoL*(usL_dot_BsL-uss_dot_Bss)*sgn(By);   ussL[4] = usR[4] + rtrhoR*(usR_dot_BsR-uss_dot_Bss)*sgn(By);

    if (SL > 0) {
        for (int i = 0; i < 8; i++)
            G[i] = GL[i];
        for (int i = 0; i < 8; i++) {
            if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
            if (fabsl(G[i]) < 1e-12) G[i] = 0;
        }
        return true;
    }
    if (SL <= 0 && SsL >= 0) {
        for (int i = 0; i < 8; i++)
            G[i] = GL[i] + SL*(usL[i]-uL[i]);
        for (int i = 0; i < 8; i++) {
            if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
            //if (fabsl(G[i]) < 1e-12) G[i] = 0;
        }
        return true;
    }
    if (SM >= 0 && SsL <= 0) {
        for (int i = 0; i < 8; i++)
            G[i] = GL[i] + SL*(usL[i]-uL[i])  +  SsL*(ussL[i]-usL[i]);
        for (int i = 0; i < 8; i++) {
            if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
            //if (fabsl(G[i]) < 1e-12) G[i] = 0;
        }
        return true;
    }

    if (SM <= 0 && SsR >= 0) {
        for (int i = 0; i < 8; i++)
            G[i] = GR[i] + SR*(usR[i]-uR[i]) + SsR*(ussR[i]-usR[i]);
        for (int i = 0; i < 8; i++) {
            if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
            //if (fabsl(G[i]) < 1e-12) G[i] = 0;
        }
        return true;
    }
    if (SR >= 0 && SsR <= 0) {
        for (int i = 0; i < 8; i++)
            G[i] = GR[i] + SR*(usR[i]-uR[i]);
        for (int i = 0; i < 8; i++) {
            if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
            //if (fabsl(G[i]) < 1e-12) G[i] = 0;
        }
        return true;
    }
    if (SR < 0) {
        for (int i = 0; i < 8; i++)
            G[i] = GR[i];
        for (int i = 0; i < 8; i++) {
            if (isnan(G[i])) { cout << "Y-HLLD: G[" << i << "] is NaN !"; return false; }
            //if (fabsl(G[i]) < 1e-12) G[i] = 0;
        }
        return true;
    }


    return false;
}




