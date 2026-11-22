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

pair<int,int> year_era(const std::string fName=""){
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
    std::string era = yr;
    yr.erase(4,6);
    era.erase(0,5);
    try{
        return std::make_pair(stoi(yr),stoi(era));
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
