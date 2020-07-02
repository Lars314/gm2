int firstSteps()
{

	//open up the file by name and list its top level contents
	TFile *f = new TFile("./gm2offline_ana_22530535_16119.00442.root"); //you can also specify options when opening a file: new, read only, write only, etc.
	f->ls();

	//change directory in the file and list again
	f->cd("farline");
	f->ls();

	//get the tree which is in the farline folder
	TTree *t = (TTree*)f->Get("farline/eventTree");
	t->Print(); //print the contents of the tree
	
	// lets get a histogram from the file
	TH1D *allCaloEnergies = (TH1D*)f->Get("farline/allCaloEnergies")->Clone("allCaloEnergies");
	allCaloEnergies->SetDirectory(0);	//set the directory to 0 so that the histogram we cloned is local
	TCanvas *c = new TCanvas(); 		//create a new canvas 
	allCaloEnergies->Draw();		//draw the histogram
	c->SetLogy();				//set the y axis to a log scale
	c->Draw();				//draw the canvas

	//create a clone of the first histogram to do our manipulations on. This will allow us to compare before/after
	TH2D* allCaloEnergies2 = (TH2D*)allCaloEnergies->Clone("allCaloEnergies2");
	double integral = allCaloEnergies->Integral(); //take the integral of the histogram. By default this encompasses all bins, but a range can be specified if you want.
	std::cout << "Integral under the curve: " << integral << std::endl; // write this out to console.
	allCaloEnergies2->Scale(1./integral);	//scale the histogram based on the result above
	std::cout << "New integral: " << allCaloEnergies2->Integral() << std::endl;
	TCanvas* c2 =  new TCanvas();
	allCaloEnergies2->Draw("hist");
	c2->SetLogy();
	c2->Draw();

	//now lets draw them together and see how they comapre
	TCanvas* c3 =  new TCanvas();
	allCaloEnergies->Draw("hist");
	allCaloEnergies->GetYaxis()->SetRangeUser(0.0001,100000000); // extend the range on the first histogram drawn so that we can see both of them
	allCaloEnergies2->SetLineColor(2);
	allCaloEnergies2->Draw("hist same");	// same means we can display the second histogram on the same axes instead of overwriting the first.
	c3->SetLogy();
	c3->Draw();


	// now lets draw the energy spectrum
	TH2D* h = new TH2D("h","Energy vs Time; Time [#mu s]; Energy [MeV]",4700,0,700,600,0,6000);
	t->Draw("energy:time*1.25/1000.>>h","","goff");
	
	//TCanvas* c4 = new TCanvas();
	TCanvas* c4 = g2Canvas2D("c",800); //this command only works with the g-2 root style installed.
	h->Draw("colz");
	c4->SetLogz();
	c4->Draw();

	return 0; // need to return something when writing C++ root macros if you make the function an int.
}
