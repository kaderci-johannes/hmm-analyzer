// LXSW - Luminosity, Cross-section, Weight

Double_t genEventSumw(TFile *f){
	TTree *runs = (TTree*)f->Get("Runs");
	Double_t out;
	runs->SetBranchAddress("genEventSumw",&out);
	runs->GetEntry(0);
	return out;
}

std::vector<std::string> runCmd(const char* cmd){
	array<char, 128> buffer;
	std::vector<std::string> result;
	FILE* pipe = popen(cmd, "r");
	if (!pipe) return {"popen failed!"};
	while (fgets(buffer.data(), buffer.size(), pipe) != nullptr) result.push_back(buffer.data());
	pclose(pipe);
	return result;
}

void LXSW(const char query[128] = "", double lumi = -1., double xsec = -1.){
	char cmd[128];
	Double_t output = 0.;
	sprintf(cmd, "ls -v ./HTCondor/lists/ | grep \"%s.*\\.list\"",query);
	cout<<endl<<"Search query : "<<cmd<<endl<<endl;
	for(string fName_list : runCmd(cmd)) cout<<fName_list;
	cout<<endl;
	for(string fName_list : runCmd(cmd)){
		if (!fName_list.empty() && fName_list.back() == '\n') fName_list.pop_back();
		string fName_root, fName_path=string("./HTCondor/lists/") + fName_list;
		ifstream f_list(fName_path);
		cout<<fName_path<<": "<<endl;
		int i=0;
		while(getline(f_list,fName_root)){
			i++;
			TFile *f = TFile::Open(fName_root.c_str(),"r");
			Double_t sumw = genEventSumw(f);
			output += sumw;
			cout<<"Sum_w_"<<i<<": "<<sumw<<endl;
			f->Close();
		}
		cout<<endl;
		f_list.close();
	}
	ofstream outfile;
	char fName_output[128];
	sprintf(fName_output,"%s%s.lxsw","./HTCondor/lists/",query);
	outfile.open(fName_output);
	cout<<"Lumi:  "<<lumi<<endl<<"xsec:  "<<xsec<<endl<<"Total: "<<output<<endl;
	outfile<<lumi<<endl<<xsec<<endl<<setprecision(25)<<output<<endl;
	outfile.close();
	return;
}

//	"root://cmsxrootd.fnal.gov//store/mc/Run3Summer22NanoAODv12/GluGluHto2Mu_M-125_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/130X_mcRun3_2022_realistic_v5-v3/50000/154d3c87-08bf-4277-ba6f-d798a21fa339.root"
