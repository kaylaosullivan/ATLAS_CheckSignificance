#include "smoother.cxx"

void testSmoother(){
    int plot = 0;
    bool absval = false;
    bool runMultiple = false;
    int numBins = 0;
    int lowLim = 0;
    int highLim = 0;

    TChain *t = new TChain("KsValidation");
    t->Add("/home/steben/projects/rrg-steven/steben/MultiquarkSearch/newtril/run/dataoct26/user.kosulliv.00267358.physics_MinBias.recon.AOD.r10170_tid13005594_00.11_hist.210140925/*.root"); 

    if (plot==0){ // tetraquark
        numBins =126; 
        lowLim = 1500;
        highLim =7000;
    }else if (plot==1){ //kaon
        numBins = 126;
        lowLim = 420;
        highLim = 580;
    }
    
    TH1D *hist1 = new TH1D ("hist1","Tetraquark (K_{s} + K_{s}) Invariant Mass", numBins, lowLim, highLim);
    
    if (plot==0){
        t->Draw("invMKK>>hist1");   
        if (runMultiple) runMultipleSigTests(hist1, 1,4,2,6,"m_{K_{s}K_{s}} [MeV]",absval);
        else significanceTest(hist1,3,3,"m_{K_{s}K_{s}} [MeV]",absval);
    }else if (plot ==1){
        t->Draw("recMass>>hist1", "RecDelR>4 && RecCosTheta>0.9998 && RecMassLambda>1125");
        if (runMultiple) runMultipleSigTests(hist1, 1,4,2,6,"m_{#pi^{+}#pi^{-}} [MeV]",absval);
        else significanceTest(hist1,2,3,"m_{#pi^{+}#pi^{-}} [MeV]",absval);
    }
   
}  
