#pragma once
#include <iostream>

#define SET_SCALAR_BRANCH(TYPE, NAME) \
    TYPE NAME{}; \
    const Bool_t has_##NAME = ((ROOT_CHAIN)->GetBranch(#NAME) != nullptr); \
    if (has_##NAME){ \
        const char* n_##TYPE = (ROOT_CHAIN)->GetLeaf(#NAME)->GetTypeName(); \
        const Bool_t is_##TYPE = (std::strcmp(n_##TYPE, #TYPE) == 0); \
        if(is_##TYPE){ \
            (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
            std::cout << "The branch " #NAME "<" << n_##TYPE << ">, is set.\n"; \
        } \
        else{ \
            std::cout << "Warning: " #TYPE " does not match with the branch " #NAME "<" << n_##TYPE << ">, branch was not set.\n"; \
        } \
    } \
    else{ \
        std::cout << "Warning: Branch \"" #NAME "\" does not exist, branch was not set.\n"; \
    } \
    static_assert(true, "")

#define SET_ARRAY_BRANCH(TYPE, NAME) \
    TYPE NAME[100]{}; \
    const Bool_t has_##NAME = ((ROOT_CHAIN)->GetBranch(#NAME) != nullptr); \
    if (has_##NAME){ \
        const char* n_##TYPE = (ROOT_CHAIN)->GetLeaf(#NAME)->GetTypeName(); \
        const Bool_t is_##TYPE = (std::strcmp(n_##TYPE, #TYPE) == 0); \
        if(is_##TYPE){ \
            (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
        } \
        else{ \
            std::cerr << "Warning: " #TYPE " does not match with the branch " #NAME "<" << n_##TYPE << ">, branch was not set.\n"; \
        } \
    } \
    else{ \
        std::cerr << "Warning: Branch \"" #NAME "\" does not exist, branch was not set.\n"; \
    } \
    static_assert(true, "")

#define SET_BRANCH_INT(NAME) \
    Int_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_UINT(NAME) \
    UInt_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_ULONG64(NAME) \
    ULong64_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_FLOAT(NAME) \
    Float_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_DOUBLE(NAME) \
    Double_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_BOOL(NAME) \
    Bool_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_UCHAR(NAME) \
    UChar_t NAME; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

    /// TODO
    /// customize further
#define SET_BRANCH_SHORT_ARRAY(NAME) \
    Short_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_USHORT_ARRAY(NAME) \
    UShort_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_INT_ARRAY(NAME) \
    Int_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_FLOAT_ARRAY(NAME) \
    Float_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_DOUBLE_ARRAY(NAME) \
    Double_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_BOOL_ARRAY(NAME) \
    Bool_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")

#define SET_BRANCH_UCHAR_ARRAY(NAME) \
    UChar_t NAME[100]; \
    (ROOT_CHAIN)->SetBranchAddress(#NAME, &NAME); \
    static_assert(true, "")
