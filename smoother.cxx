#include <math.h>
#include <iostream>

using namespace std;

bool rounding = false;
bool absval = false;

// Function to smooth data numSmooths amount of times
// Each point is smoothed by: x_inew = (x_(i-1) + 2*x_i + x_(i+1))/4
// First and last points are smoothed by: x_1new = (2*x_1 + x_2)/3  and x_Nnew = (x_(N-1) + 2*x_N)/3
TH1D* smooth(TH1D* data, int numSmooths){
    
    TH1D *smoothed_data = (TH1D*)data->Clone(); // Clone data
    smoothed_data->SetName("smoothed_data");
    int x;

    int smoothedSize = smoothed_data->GetSize()-2; // subtract 2 to remove under/overflow bins
    
    // loop through each element and calculate new smoothed point
    for(int i = 1; i <= numSmooths; i++){
	for (int j = i+1; j <= smoothedSize-i; j++){
	    
	    if (j==1) // first point
		    x = (2*smoothed_data->GetBinContent(j) + smoothed_data->GetBinContent(j+1))/3; 
	    else if (j == smoothedSize) // second point
		    x = (smoothed_data->GetBinContent(j-1) + 2*smoothed_data->GetBinContent(j))/3;
	    else  // other points
 	        x = (smoothed_data->GetBinContent(j-1) + 2*smoothed_data->GetBinContent(j) + smoothed_data->GetBinContent(j+1))/4;
        
        smoothed_data->SetBinContent(j,x);
        }
    }
    return smoothed_data;
}

// Calculates the significance at each point by: (x_old - x_smoothed)/std, where std = sqrt(x_old)
// Note that there cannot be a bin with zero content
TH1D* significance(TH1D* data, TH1D* smoothed_data, int mergeBins){
    
    TH1D *sig = (TH1D*)smoothed_data->Clone();
    sig->SetName("sig");
    
    // Initialize values
    int smoothedSize = smoothed_data->GetSize()-2;
    double std;
    double diff;
    double signif;    
    
    for(int i = 1; i <= smoothedSize; i++){
        std = sqrt(data->GetBinContent(i)); // Calculate standard deviation 
        
    	if (std==0){ // exit if there is zero bin content
    	     cout<<("There exists a bin with 0 content. Please rebin original histogram.")<< endl;
    	     exit(0);    
    	}	
           
        // Calculate the difference and significance
        diff = data->GetBinContent(i) - smoothed_data->GetBinContent(i);
        sig->SetBinContent(i,diff/std);
    }
    
    // Merge specified number of bins
    sig->Rebin(mergeBins);
 
    // Take absolute value or round if specified 
    if (rounding == true || absval == true){
	    for (int i = 1; i<=smoothedSize; i++){
      	     if (rounding == true) sig->SetBinContent(i,round(sig->GetBinContent(i)));
    	     if (absval == true) sig->SetBinContent(i,fabs(sig->GetBinContent(i)));
        }
    }
       
    return sig;
}

// Test for significance: plots smoothed data and significance on separate canvases
// option for printInfo puts number of smooths and merged bins on plot
// option for drawLine plots curve instead of points for smoothed data
void significanceTest(TH1D* hist, int numSmooths, int mergeBins, const char* title, bool absVal=false,  bool drawLine=false, bool roundVal=false, bool printInfo=true){
    rounding = roundVal;
    absval = absVal;     

    TH1D *smoothed = (TH1D*)smooth(hist, numSmooths);
    TH1D *signif = (TH1D*)significance(hist, smoothed,mergeBins);   

    // Create canvases
    TCanvas *c1 = new TCanvas("c1"); 
    TCanvas *c2 = new TCanvas("c2");    

    // Histogram settings
    c1->cd();
    hist->SetFillColor(kCyan-10);
    signif->SetFillColor(kMagenta-10);
    smoothed->SetMarkerColor(kRed+1);;    
    hist->SetLineColor(kBlue);
    smoothed->SetLineColor(kRed);
    smoothed->SetLineWidth(2);
    smoothed->SetMarkerSize(1);

    // Histogram titles
    hist->GetXaxis()->SetTitle(title);
    hist->GetYaxis()->SetTitle("Number of Events per Bin");
    signif->GetYaxis()->SetTitle("Standard Deviations from Smoothed Data");
    signif->GetXaxis()->SetTitle(title);    
    signif->GetYaxis()->SetTitleSize(.038);  

    // Draw histogram
    hist->Draw();
    smoothed->Draw("PSAME"); // Add "C" for smooth curve
    if (drawLine == true) smoothed->Draw("CSAME");

    // Add Legend
    TLegend *leg = new TLegend(.65,.78,.93,.9);
    leg->AddEntry(hist, "Data","f");   // p = show marker
    if (drawLine == true) leg->AddEntry(smoothed,("Smoothed Data ("+std::to_string(numSmooths)+"x)").c_str(),"lp");
    else leg->AddEntry(smoothed, ("Smoothed Data ("+std::to_string(numSmooths)+"x)").c_str(), "p");          // f = show fill 
    leg->Draw("same");
    c1->Update();
 
    // Draw Significance
    c2->cd();
    signif->Draw("HIST");

    // Print information on plot
    if (printInfo == true){
    	TPaveText *pt = new TPaveText(0.8,0.9,0.99,0.99,"NDC");
      	pt->SetFillColor(0);
    	pt->SetBorderSize(1);
    	pt->AddText(("Data Smoothed "+std::to_string(numSmooths)+"x").c_str());
    	pt->AddText((std::to_string(mergeBins)+" Bins Merged").c_str());
    	pt->Draw("same");
    }
}

// Runs multiple significance tests for varying number of smooths and merged bins
void runMultipleSigTests(TH1D* hist, int lowSmooths, int highSmooths, int lowBins, int highBins, const char* title, bool absVal=false, bool round=false){
    
    TPaveText *pt = new TPaveText(0,0,1,1,"NDC");
    pt->SetFillColor(0);

    // Initialize values
    int x = highSmooths - lowSmooths + 1;
    int y = highBins - lowBins + 1;
    int numj = 1;
    int numi = 1;
    int yfixed = y;
    int numBins = hist->GetSize() -2;
    bool isSig = false;

    // Rows of plots excluding indivisible numbers of smoothed bins
    for (int k=lowBins; k <= highBins;k++ ) {
	    if (numBins%k != 0) yfixed--;
    }
   
   // Create and divide canvases
   // For smoothed data
   TCanvas *c1 = new TCanvas("c1","",1080,700);
   c1->Divide(ceil(sqrt(x)),ceil(x/ceil(sqrt(x))),0.001,0.001);

    // For significance plots
    TCanvas *c2 = new TCanvas("c2","",1080,700);
    c2->Divide(x,yfixed,0.001,0.001);

    // Make plots for each combination of merged bins and smooths
    for(int j=lowBins; j<=highBins;j++){
        for(int i=lowSmooths;i<=highSmooths;i++){
    	    
    	    // Skip if number of bins is indivisible by number of bins to merge
    	    if (numBins%j != 0){
        		cout << "Skipping bin merge by " << j << " since total number of bins is indivisible by " << j << endl;
        		break;
	        }
	        
	        TH1D *smoothed = (TH1D*)smooth(hist, i);
	        TH1D *signif = (TH1D*)significance(hist, smoothed, j);
            
            // Print significance info (observation for sigma >=3, discovery for sigma >=5)
		    double maxVal = signif->GetBinContent(signif->GetMaximumBin()); 
    		if (maxVal >= 5){
    			cout<<"Possible discovery when smoothed "<<i<<" times, with "<<j<<" bins merged."<<endl;
    			isSig = true;
    	  	}else if (maxVal >= 3){
    			cout<<"Possible observation when smoothed "<<i<<" times, with "<<j<<" bins merged."<<endl;
    			isSig = true;
    		}
            
            // Histogram Settings
    	    hist->SetFillColor(kCyan-10);
            signif->SetFillColor(kMagenta-10);
        	smoothed->SetMarkerColor(kRed+1);
        	smoothed->SetLineColor(kRed);
    	    smoothed->SetMarkerSize(1);
    
            // Histogram titles
            hist->GetXaxis()->SetTitle(title);
            hist->GetYaxis()->SetTitle("Number of Events per Bin");
        	signif->GetYaxis()->SetTitle("Standard Deviations from Smoothed Data");
        	signif->GetXaxis()->SetTitle(title);
        	signif->GetYaxis()->SetTitleSize(.038);
    
            // Draw smoothed data only when it changes
            if (j == lowBins){
    	        c1->cd(numj);
        	    hist->Draw();
        	    smoothed->Draw("PSAME"); // Add "C" for smooth curve
            	            	    
            	// Add legend
            	TLegend *leg = new TLegend(.6,.8,.93,.9);
        	    leg->AddEntry(hist, "Data","f");   // p = show marker
           	    leg->AddEntry(smoothed, ("Smoothed Data ("+std::to_string(i)+"x)").c_str(), "p");     // f = show fill
          	    leg->Draw("same");
            	c1->Update();
        		numj++;
            }

            // Draw significance histogram
	        c2->cd(numi);
            signif->Draw("HIST");
            
            // Draw smoothing and merging info
            TPaveText *pt = new TPaveText(0.6,0.8,0.99,0.99,"NDC");
            pt->SetFillColor(0);
            pt->SetBorderSize(1);
            pt->AddText(("Data Smoothed "+std::to_string(i)+"x").c_str());
            pt->AddText((std::to_string(j)+" Bins Merged").c_str());
            pt->Draw("same");

	        numi++;
	    }	
    }
    // If no observations or disoveries found
    if (!isSig) cout<<"Nothing significant was found in the data." << endl;
}
