void int_assign(bool isSigned, Int_t& nItem, Int_t& nItemSigned, UInt_t& nItemUnsigned, bool isMC=true){
    if(!isMC){
        nItem=0;
        return;
    }
    if(isSigned){
        nItem=nItemSigned;
    }
    else{
        nItem=nItemUnsigned;
    }
}

void int_assign(bool isSigned, Int_t& nItem, Int_t& nItemSigned, UChar_t& nItemUnsigned, bool isMC=true){
    if(!isMC){
        nItem=0;
        return;
    }
    if(isSigned){
        nItem=nItemSigned;
    }
    else{
        nItem=nItemUnsigned;
    }
}

pair<int,int> year_semester(const std::string fName=""){
    std::string yr = fName;
    size_t posmc = yr.find("mc_");
    size_t posdata = yr.find("data_");
    size_t poseff = yr.find("eff_");
    if (posmc != std::string::npos){
        yr.erase(0,posmc+3);
    }
    else if(posdata != std::string::npos){
        yr.erase(0,posdata+5);
    }
    else if(poseff != std::string::npos){
        yr.erase(0,poseff+4);
    }
    yr.erase(6,-1);
    std::string semester = yr;
    yr.erase(4,6);
    semester.erase(0,5);
    try{
        return std::make_pair(stoi(yr),stoi(semester));
    }
    catch(std::invalid_argument &e1){
        try{
            return std::make_pair(stoi(yr),1);
        }
        catch(std::invalid_argument &e2){
            throw std::runtime_error("Invalid year on --input filename.");
        }
    }
}

std::string era(const std::string fName=""){
    std::string s = fName;
    size_t posMuon0 = s.find("Muon0_");
    size_t posMuon1 = s.find("Muon1_");
    size_t posMuon = s.find("Muon_");
    if(posMuon0 != std::string::npos){
        s.erase(0,posMuon0+6);
        s.erase(1,-1);
    }
    else if(posMuon1 != std::string::npos){
        s.erase(0,posMuon1+6);
        s.erase(1,-1);
    }
    else if(posMuon != std::string::npos){
        s.erase(0,posMuon+5);
        s.erase(1,-1);
    }
    else{
        throw std::runtime_error("era() - Invalid --input filename for ntuple mode.");
    }
    return s;
}

std::string capture_string(const std::string fName="", const std::string search="", const int n=1){
    std::string s = fName;
    size_t posMuon = s.find(search);
    if(posMuon != std::string::npos){
        s.erase(0,posMuon+search.size());
        s.erase(n,-1);
    }
    else{
        throw std::runtime_error("capture_string() - Invalid --input filename.");
    }
    return s;
}

std::string sample_id(const std::string fName=""){
    std::string s = fName;
    size_t posMuon = s.find("Muon");
    size_t posggH = s.find("ggH");
    if(posMuon != std::string::npos){
        s.erase(0,posMuon);
        s.erase(5,-1);
    }
    else if(posggH != std::string::npos){
        s.erase(0,posggH+4);
        s.erase(3,-1);
    }
    else{
        throw std::runtime_error("Invalid --input filename for ntuple mode.");
    }
    return s;
}
