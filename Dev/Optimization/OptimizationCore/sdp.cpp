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
        n_ += C.blocks[index].blocksize;
        ++index;
    }
    C_ =(void*)Cptr;
}

void SDP::setAs(const std::vector<std::vector<arma::mat>>& As)
{
    k_ = As.size();
    struct constraintmatrix* constraints = (struct constraintmatrix *)malloc((k_+1)*sizeof(struct constraintmatrix));
    struct sparseblock* blockptr;
    if(NULL==constraints)
    {
        std::cerr<<"Couldn't allocate storage for constraints"<<std::endl;
    }
    int index = 1;
    std::vector<std::vector<arma::mat>>::const_iterator asiter;

    for(asiter=As.cbegin();asiter!=As.cend();++asiter)
    {
        constraints[index].blocks=NULL;
        //the link list can start with block at tail
        int bindex = asiter->size();
        std::vector<arma::mat>::const_reverse_iterator aiter;
        for(aiter=asiter->crbegin();aiter!=asiter->crend();++aiter)
        {
            if(!aiter->empty())
            {
                blockptr=(struct sparseblock*)malloc(sizeof(struct sparseblock));
                if(NULL==blockptr)
                {
                    assert(blockptr);
                }
                blockptr->blocknum = bindex;
                blockptr->blocksize = aiter->n_rows;
                blockptr->constraintnum = index;
                blockptr->next = NULL;
                blockptr->nextbyblock=NULL;
                blockptr->numentries=0;
                if(aiter->n_rows==aiter->n_cols)
                {
                    int r,c;
                    for(r=0;r<aiter->n_rows;++r)
                        for(c=r;c<aiter->n_cols;++c)
                        {
                           if((*aiter)(r,c))
                           {
                               blockptr->numentries++;
                           }
                        }
                }else{//diag block
                    for(int r=0;r<aiter->n_rows;++r)
                    {
                        if((*aiter)(r))
                        {
                            blockptr->numentries++;
                        }
                    }
                }
                blockptr->entries=(double*)malloc((blockptr->numentries+1)*sizeof(double));
                if(NULL==blockptr->entries)
                {
                    assert(blockptr->entries);
                }
                blockptr->iindices=(unsigned short*)malloc((blockptr->numentries+1)*sizeof(unsigned short));
                if(NULL==blockptr->iindices)
                {
                    assert(blockptr->iindices);
                }
                blockptr->jindices=(unsigned short*)malloc((blockptr->numentries+1)*sizeof(unsigned short));
                if(NULL==blockptr->jindices)
                {
                    assert(blockptr->jindices);
                }
                if(aiter->n_rows==aiter->n_cols)
                {
                    int r,c;
                    int i=1;
                    for(r=0;r<aiter->n_rows;++r)
                        for(c=r;c<aiter->n_cols;++c)
                        {
                           if((*aiter)(r,c))
                           {
                               blockptr->entries[i]=(*aiter)(r,c);
                               blockptr->iindices[i]=r+1;
                               blockptr->jindices[i]=c+1;
                               ++i;
                           }
                        }
                }else{//diag block
                    int i=1;
                    for(int r=0;r<aiter->n_rows;++r)
                    {
                        if((*aiter)(r))
                        {
                            blockptr->entries[i]=(*aiter)(r);
                            blockptr->iindices[i]=r+1;
                            blockptr->jindices[i]=r+1;
                            ++i;
                        }
                    }
                }
                blockptr->next=constraints[index].blocks;
                constraints[index].blocks=blockptr;
            }
            --bindex;
        }
        ++index;
    }
    constraints_ = (void*)constraints;
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
    struct constraintmatrix* constraints = (struct constraintmatrix*)constraints_;
    if(NULL==C)return false;
    if(NULL==constraints)return false;
    if(NULL==b_)return false;
    X_ = (void*)(new struct blockmatrix);
    struct blockmatrix* X = (struct blockmatrix*)X_;
    if(NULL==X)return false;
    Z_ = (void*)(new struct blockmatrix);
    struct blockmatrix* Z = (struct blockmatrix*)Z_;
    if(NULL==Z)return false;
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
    arma::vec tmp(y_,k_+1,false,true);
    y = tmp.tail(k_);
}
void SDP::getX(arma::sp_mat& Xmat)
{
    if(!X_)return;
    struct blockmatrix* X = (struct blockmatrix*)X_;
    arma::uword N = 0;
    for(int bn=1;bn<=X->nblocks;++bn)
    {
        struct blockrec& block = X->blocks[bn];
        N += block.blocksize*block.blocksize;
    }
    std::vector<double> value;
    value.reserve(N);
    std::vector<arma::uword> rows;
    rows.reserve(N);
    std::vector<arma::uword> cols;
    cols.reserve(N);
    arma::uword global_r=0,global_c=0;
    for(int bn=1;bn<=X->nblocks;++bn)
    {
        struct blockrec& block = X->blocks[bn];
        if(block.blockcategory==MATRIX)
        {
            for(int dr=1;dr<=block.blocksize;++dr)
                for(int dc=1;dc<=block.blocksize;++dc)
                {
                    if(0.0!=block.data.mat[ijtok(dr,dc,block.blocksize)])
                    {
                        value.push_back(block.data.mat[ijtok(dr,dc,block.blocksize)]);
                        rows.push_back(global_r+dr-1);
                        cols.push_back(global_c+dc-1);
                    }
                }
        }
        if(block.blockcategory==DIAG)
        {
            for(int di=1;di<=block.blocksize;++di)
            {
                if(block.data.vec[di]!=0.0)
                {
                    value.push_back(block.data.vec[di]);
                    rows.push_back(global_r+di-1);
                    cols.push_back(global_c+di-1);
                }
            }
        }
        global_r += block.blocksize;
        global_c += block.blocksize;
    }
    arma::umat location(2,rows.size());
    location.row(0) = arma::urowvec(rows.data(),rows.size(),false,true);
    location.row(1) = arma::urowvec(cols.data(),rows.size(),false,true);
    Xmat = arma::sp_mat(location,arma::vec(value.data(),value.size(),false,true));
}
void SDP::getZ(arma::sp_mat& Zmat)
{
    if(!Z_)return;
    struct blockmatrix* Z = (struct blockmatrix*)Z_;
    arma::uword N = 0;
    for(int bn=1;bn<=Z->nblocks;++bn)
    {
        struct blockrec& block = Z->blocks[bn];
        N += block.blocksize*block.blocksize;
    }
    std::vector<double> value;
    value.reserve(N);
    std::vector<arma::uword> rows;
    rows.reserve(N);
    std::vector<arma::uword> cols;
    cols.reserve(N);
    arma::uword global_r=0,global_c=0;
    for(int bn=1;bn<=Z->nblocks;++bn)
    {
        struct blockrec& block = Z->blocks[bn];
        if(block.blockcategory==MATRIX)
        {
            for(int dr=1;dr<=block.blocksize;++dr)
                for(int dc=1;dc<=block.blocksize;++dc)
                {
                    if(0.0!=block.data.mat[ijtok(dr,dc,block.blocksize)])
                    {
                        value.push_back(block.data.mat[ijtok(dr,dc,block.blocksize)]);
                        rows.push_back(global_r+dr-1);
                        cols.push_back(global_c+dc-1);
                    }
                }
        }
        if(block.blockcategory==DIAG)
        {
            for(int di=1;di<=block.blocksize;++di)
            {
                if(block.data.vec[di]!=0.0)
                {
                    value.push_back(block.data.vec[di]);
                    rows.push_back(global_r+di-1);
                    cols.push_back(global_c+di-1);
                }
            }
        }
        global_r += block.blocksize;
        global_c += block.blocksize;
    }
    arma::umat location(2,rows.size());
    location.row(0) = arma::urowvec(rows.data(),rows.size(),false,true);
    location.row(1) = arma::urowvec(cols.data(),rows.size(),false,true);
    Zmat = arma::sp_mat(location,arma::vec(value.data(),value.size(),false,true));
}

void SDP::debug_prob(char* path)
{
    struct blockmatrix* Cptr = (struct blockmatrix*)C_;
    struct blockmatrix& C = *Cptr;
    struct constraintmatrix* constraints = (struct constraintmatrix*)constraints_;
    write_prob(path,7,2,C,b_,constraints);
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
                if(ptr->entries)free(ptr->entries);
                if(ptr->iindices)free(ptr->iindices);
                if(ptr->jindices)free(ptr->jindices);
                oldptr=ptr;
                ptr=ptr->next;
                free(oldptr);
            };
        };
        free(constraints);
    };
}
}

