#include "sdp.h"
#include <strstream>
#include <csdp/declarations.h>
#include <assert.h>
namespace  Optimization {
SDP::SDP():
    b_(NULL),
    constraints_(NULL),
    y_(NULL),
    C_(NULL),
    X_(NULL),
    Z_(NULL)
{
    ;
}

void SDP::setC(const std::vector<arma::mat>&iC)
{
    n_ = 0;
    struct blockmatrix* Cptr = new struct blockmatrix;
    struct blockmatrix& C = *Cptr;
    C.nblocks = 0;
    C.blocks=(struct blockrec*)malloc((iC.size()+1)*sizeof(struct blockrec));
    if (C.blocks == NULL)
    {
        std::cerr<<"Couldn't allocate storage for C!"<<std::endl;
        delete Cptr;
        return;
    };
    //set up blocks
    std::vector<arma::mat>::const_iterator citer;
    int index = 1;
    for(citer=iC.cbegin();citer!=iC.cend();++citer)
    {
        const arma::mat& iblock = *citer;
        C.blocks[index].blocksize = iblock.n_rows;
        int size = C.blocks[index].blocksize;
        if(iblock.n_rows==iblock.n_cols)
        {
            C.blocks[index].blockcategory = MATRIX;
            C.blocks[index].data.mat=(double*)malloc(
                        size*size*
                        sizeof(double)
                        );
            if (C.blocks[index].data.mat == NULL)
            {
                std::cerr<<"Couldn't allocate storage for "
                         <<"C.block["<<index<<"] with size:"
                         <<C.blocks[index].blocksize<<" and type: MATRIX"<<std::endl;
                free_mat(C);
                delete Cptr;
                n_ = 0;
                return;
            };
            arma::mat tmp(C.blocks[index].data.mat,size,size,false,true);
            tmp=iblock;
        }else{
            assert(iblock.n_cols==1);
            C.blocks[index].blockcategory = DIAG;
            C.blocks[index].data.vec=(double *)malloc((size+1)*sizeof(double));
            if (C.blocks[index].data.vec == NULL)
            {
                std::cerr<<"Couldn't allocate storage for "
                         <<"C.block["<<index<<"] with size:"
                         <<C.blocks[index].blocksize<<" and type: DIAG"<<std::endl;
                free_mat(C);
                delete Cptr;
                n_ = 0;
                return;
            };
            arma::vec tmp(C.blocks[index].data.vec,size+1,false,true);
            tmp.tail(size) = iblock;
        }
        ++C.nblocks; //increase nblocks only when the allocate is successfully done
        ++index;
    }
    C_ =(void*)Cptr;
}

void SDP::setAs(const std::vector<std::vector<arma::mat>>& As)
{
    k_ = As.size();
    ;
}

void SDP::setb(const arma::vec& b)
{
    b_=(double *)malloc((b.size()+1)*sizeof(double));
    if (b_==NULL)
    {
        std::cerr<<"Couldn't allocate storage for b with size:"<<b.size()<<std::endl;
        return;
    };
    arma::vec tmp(b_,b.size()+1,false,true);
    tmp.tail(b.size()) = b;
}

bool SDP::init()
{
    struct blockmatrix* C = (struct blockmatrix*)C_;
    struct blockmatrix* X = (struct blockmatrix*)X_;
    struct blockmatrix* Z = (struct blockmatrix*)Z_;
    struct constraintmatrix* constraints = (struct constraintmatrix*)constraints_;
    if(NULL==C)return false;
    if(NULL==constraints)return false;
    if(NULL==b_)return false;
    initsoln(n_,k_,*C,b_,constraints,X,&y_,Z);
    return true;
}

bool SDP::solve()
{
    struct blockmatrix* C = (struct blockmatrix*)C_;
    struct blockmatrix* X = (struct blockmatrix*)X_;
    struct blockmatrix* Z = (struct blockmatrix*)Z_;
    struct constraintmatrix* constraints = (struct constraintmatrix*)constraints_;
    ret_ = easy_sdp(n_,k_,*C,b_,constraints,0.0,X,&y_,Z,&pobj_,&dobj_);
    if(ret_!=0)return false;
    return true;
}

std::string SDP::info()const
{
    std::stringstream stream;
    if (ret_ == 0)
    {
        stream<<"The objective value is "<<(dobj_+pobj_)/2<<"\n";
    }
    else
    {
        stream<<"SDP failed\n";
    }
    return stream.str();
}
void SDP::gety(arma::vec& y)
{
    ;
}
void SDP::getX(arma::sp_mat& X)
{
    ;
}
void SDP::getZ(arma::sp_mat& Z)
{
    ;
}
SDP::~SDP()
{
    struct blockmatrix* C = (struct blockmatrix*)C_;
    struct blockmatrix* X = (struct blockmatrix*)X_;
    struct blockmatrix* Z = (struct blockmatrix*)Z_;
    struct constraintmatrix* constraints = (struct constraintmatrix*)constraints_;
    if(y_)free(y_);
    if(b_)free(b_);
    if(C){
        free_mat(*C);
        delete C;
        C_ = NULL;
    }
    if(X){
        free_mat(*X);
        delete X;
        X_ = NULL;
    }
    if(Z){
        free_mat(*Z);
        delete Z;
        Z_ = NULL;
    }
    int i;
    struct sparseblock *ptr;
    struct sparseblock *oldptr;
    if(constraints)
    {
        for (i=1; i<=k_; i++)
        {
            ptr=constraints[i].blocks;
            while (ptr != NULL)
            {
                free(ptr->entries);
                free(ptr->iindices);
                free(ptr->jindices);
                oldptr=ptr;
                ptr=ptr->next;
                free(oldptr);
            };
        };
        free(constraints);
    };
}
}

