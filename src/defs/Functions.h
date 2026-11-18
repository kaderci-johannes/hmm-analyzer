void int_assign(bool isSigned, Int_t& nItem, Int_t& nItemSigned, UInt_t& nItemUnsigned){
    if(isSigned){
        nItem=nItemSigned;
    }
    else{
        nItem=nItemUnsigned;
    }
}
