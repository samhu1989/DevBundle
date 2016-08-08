/*
    Copyright (c) 2013, Philipp Krähenbühl
    All rights reserved.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are met:
        * Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.
        * Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.
        * Neither the name of the Stanford University nor the
        names of its contributors may be used to endorse or promote products
        derived from this software without specific prior written permission.

    THIS SOFTWARE IS PROVIDED BY Philipp Krähenbühl ''AS IS'' AND ANY
    EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
    WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
    DISCLAIMED. IN NO EVENT SHALL Philipp Krähenbühl BE LIABLE FOR ANY
    DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
    (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
    LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
    ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "densecrf.h"
#include "CRF/permutohedral.h"
#include "CRF/util.h"
#include "CRF/pairwise.h"
#include <cmath>
#include <cstring>
#include <iostream>
#include <QRgb>
/////////////////////////////
/////  Alloc / Dealloc  /////
/////////////////////////////
DenseCRF::DenseCRF(int N, int M) : N_(N), M_(M), unary_(0) {
}
DenseCRF::~DenseCRF() {
	if (unary_)
		delete unary_;
	for( unsigned int i=0; i<pairwise_.size(); i++ )
		delete pairwise_[i];
}
DenseCRF2D::DenseCRF2D(int W, int H, int M) : DenseCRF(W*H,M), W_(W), H_(H) {
}
DenseCRF2D::~DenseCRF2D() {
}
/////////////////////////////////
/////  Pairwise Potentials  /////
/////////////////////////////////
void DenseCRF::addPairwiseEnergy (const arma::mat & features, LabelCompatibility * function, KernelType kernel_type, NormalizationType normalization_type) {
    assert( features.n_cols == N_ );
	addPairwiseEnergy( new PairwisePotential( features, function, kernel_type, normalization_type ) );
}
void DenseCRF::addPairwiseEnergy ( PairwisePotential* potential ){
	pairwise_.push_back( potential );
}
void DenseCRF2D::addPairwiseGaussian ( float sx, float sy, LabelCompatibility * function, KernelType kernel_type, NormalizationType normalization_type ) {
    arma::mat feature( 2, N_ );
	for( int j=0; j<H_; j++ )
		for( int i=0; i<W_; i++ ){
			feature(0,j*W_+i) = i / sx;
			feature(1,j*W_+i) = j / sy;
		}
	addPairwiseEnergy( feature, function, kernel_type, normalization_type );
}
void DenseCRF2D::addPairwiseBilateral ( float sx, float sy, float sr, float sg, float sb, const unsigned char* im, LabelCompatibility * function, KernelType kernel_type, NormalizationType normalization_type ) {
    arma::mat feature( 5, N_ );
	for( int j=0; j<H_; j++ )
		for( int i=0; i<W_; i++ ){
			feature(0,j*W_+i) = i / sx;
			feature(1,j*W_+i) = j / sy;
			feature(2,j*W_+i) = im[(i+j*W_)*3+0] / sr;
			feature(3,j*W_+i) = im[(i+j*W_)*3+1] / sg;
			feature(4,j*W_+i) = im[(i+j*W_)*3+2] / sb;
		}
	addPairwiseEnergy( feature, function, kernel_type, normalization_type );
}
//////////////////////////////
/////  Unary Potentials  /////
//////////////////////////////
void DenseCRF::setUnaryEnergy ( UnaryEnergy * unary ) {
	if( unary_ ) delete unary_;
	unary_ = unary;
}
void DenseCRF::setUnaryEnergy( const arma::mat & unary ) {
	setUnaryEnergy( new ConstUnaryEnergy( unary ) );
}
void  DenseCRF::setUnaryEnergy( const arma::mat & L, const arma::mat & f ) {
	setUnaryEnergy( new LogisticUnaryEnergy( L, f ) );
}
///////////////////////
/////  Inference  /////
///////////////////////
void expAndNormalize ( arma::mat & out, const arma::mat & in ) {
    out = arma::mat( in.n_rows, in.n_cols , arma::fill::zeros );
    #pragma omp parallel for
    for( int i=0; i<out.n_cols; i++ ){
        out.col(i) = in.col(i);
        out.col(i) -= arma::max(out.col(i));
        out.col(i) = arma::exp(out.col(i));
        out.col(i) /= arma::accu(out.col(i));
	}
}

void sumAndNormalize( arma::mat & out, const arma::mat & in, const arma::mat & Q ) {
    out = arma::mat( in.n_rows, in.n_cols , arma::fill::zeros );
    #pragma omp parallel for
    for( int i=0; i<in.n_cols; i++ ){
        out.col(i) = Q.col(i);
        out.col(i) *= arma::accu(in.col(i));
        out.col(i) -= in.col(i);
	}
}

arma::mat DenseCRF::inference ( int n_iterations ) const {
    arma::mat Q( M_, N_ ), tmp1, unary( M_, N_, arma::fill::zeros ), tmp2;
	if( unary_ )
		unary = unary_->get();
	expAndNormalize( Q, -unary );
	
	for( int it=0; it<n_iterations; it++ ) {
		tmp1 = -unary;
		for( unsigned int k=0; k<pairwise_.size(); k++ ) {
			pairwise_[k]->apply( tmp2, Q );
			tmp1 -= tmp2;
		}
		expAndNormalize( Q, tmp1 );
	}
	return Q;
}

arma::uvec DenseCRF::map ( int n_iterations ) const {
	// Run inference
    arma::mat Q = inference( n_iterations );
	// Find the map
	return currentMap( Q );
}
///////////////////
/////  Debug  /////
///////////////////
arma::vec DenseCRF::unaryEnergy(const arma::uvec & l) {
    assert( l.n_rows == N_ );
    arma::vec r( N_ ,arma::fill::zeros );
	if( unary_ ) {
        arma::mat unary = unary_->get();
		for( int i=0; i<N_; i++ )
			if ( 0 <= l[i] && l[i] < M_ )
				r[i] = unary( l[i], i );
	}
	return r;
}

arma::vec DenseCRF::pairwiseEnergy(const arma::uvec & l, int term) {
    assert( l.n_rows == N_ );
    arma::vec r( N_ ,arma::fill::zeros );
	if( term == -1 ) {
		for( unsigned int i=0; i<pairwise_.size(); i++ )
			r += pairwiseEnergy( l, i );
		return r;
	}
    arma::mat Q( M_, N_ );
	// Build the current belief [binary assignment]
	for( int i=0; i<N_; i++ )
		for( int j=0; j<M_; j++ )
			Q(j,i) = (l[i] == j);
	pairwise_[ term ]->apply( Q, Q );
	for( int i=0; i<N_; i++ )
		if ( 0 <= l[i] && l[i] < M_ )
			r[i] =-0.5*Q(l[i],i );
		else
			r[i] = 0;
	return r;
}

arma::mat DenseCRF::startInference() const{
    arma::mat Q( M_, N_, arma::fill::zeros );
	// Initialize using the unary energies
	if( unary_ )
		expAndNormalize( Q, -unary_->get() );
	return Q;
}
void DenseCRF::stepInference( arma::mat & Q, arma::mat & tmp1, arma::mat & tmp2 ) const{
    tmp1 = arma::mat( Q.n_rows, Q.n_cols , arma::fill::zeros );
	if( unary_ )
		tmp1 -= unary_->get();
	
	// Add up all pairwise potentials
	for( unsigned int k=0; k<pairwise_.size(); k++ ) {
		pairwise_[k]->apply( tmp2, Q );
		tmp1 -= tmp2;
	}
	
	// Exponentiate and normalize
	expAndNormalize( Q, tmp1 );
}
arma::uvec DenseCRF::currentMap( const arma::mat & Q ) const{
    arma::uvec r(Q.n_cols);
	// Find the map
    #pragma omp parallel for
	for( int i=0; i<N_; i++ ){
        arma::uword m;
        Q.col(i).max(m);
		r[i] = m;
	}
	return r;
}

// Compute the KL-divergence of a set of marginals
double DenseCRF::klDivergence( const arma::mat & Q ) const {
	double kl = 0;
	// Add the entropy term
    for( int i=0; i<Q.n_cols; i++ )
        for( int l=0; l<Q.n_rows; l++ )
            kl += Q(l,i)*std::log(std::max( Q(l,i), double(1e-20f)) );
	// Add the unary term
	if( unary_ ) {
        arma::mat unary = unary_->get();
        for( int i=0; i<Q.n_cols; i++ )
            for( int l=0; l<Q.n_rows; l++ )
				kl += unary(l,i)*Q(l,i);
	}
	
	// Add all pairwise terms
    arma::mat tmp;
	for( unsigned int k=0; k<pairwise_.size(); k++ ) {
		pairwise_[k]->apply( tmp, Q );
        kl += arma::accu(Q%tmp);
	}
	return kl;
}

// Gradient computations
double DenseCRF::gradient( int n_iterations, const ObjectiveFunction & objective, arma::vec * unary_grad, arma::vec * lbl_cmp_grad, arma::vec * kernel_grad) const {
	// Run inference
    std::vector<arma::mat> Q(n_iterations+1);
    arma::mat tmp1, unary( M_, N_ ), tmp2;
	unary.fill(0);
	if( unary_ )
		unary = unary_->get();
	expAndNormalize( Q[0], -unary );
	for( int it=0; it<n_iterations; it++ ) {
		tmp1 = -unary;
		for( unsigned int k=0; k<pairwise_.size(); k++ ) {
			pairwise_[k]->apply( tmp2, Q[it] );
			tmp1 -= tmp2;
		}
		expAndNormalize( Q[it+1], tmp1 );
	}
	
	// Compute the objective value
    arma::mat b( M_, N_ );
	double r = objective.evaluate( b, Q[n_iterations] );
	sumAndNormalize( b, b, Q[n_iterations] );

	// Compute the gradient
	if(unary_grad && unary_)
		*unary_grad = unary_->gradient( b );
	if( lbl_cmp_grad )
		*lbl_cmp_grad = 0*labelCompatibilityParameters();
	if( kernel_grad )
		*kernel_grad = 0*kernelParameters();
	
	for( int it=n_iterations-1; it>=0; it-- ) {
		// Do the inverse message passing
		tmp1.fill(0);
		int ip = 0, ik = 0;
		// Add up all pairwise potentials
		for( unsigned int k=0; k<pairwise_.size(); k++ ) {
			// Compute the pairwise gradient expression
			if( lbl_cmp_grad ) {
                arma::mat pg = pairwise_[k]->gradient( b, Q[it] );
                lbl_cmp_grad->subvec( ip, ip + pg.n_rows - 1 ) += pg;
                ip += pg.n_rows;
			}
			// Compute the kernel gradient expression
			if( kernel_grad ) {
                arma::mat pg = pairwise_[k]->kernelGradient( b, Q[it] );
                kernel_grad->subvec( ik, ik + pg.n_rows - 1 ) += pg;
                ik += pg.n_rows;
			}
			// Compute the new b
			pairwise_[k]->applyTranspose( tmp2, b );
			tmp1 += tmp2;
		}
        sumAndNormalize( b, tmp1%Q[it], Q[it] );
		
		// Add the gradient
		if(unary_grad && unary_)
			*unary_grad += unary_->gradient( b );
	}
	return r;
}

arma::vec DenseCRF::unaryParameters() const {
	if( unary_ )
		return unary_->parameters();
    return arma::vec();
}

void DenseCRF::setUnaryParameters( const arma::vec & v ) {
	if( unary_ )
		unary_->setParameters( v );
}

arma::vec DenseCRF::labelCompatibilityParameters() const {
    std::vector< arma::vec > terms;
	for( unsigned int k=0; k<pairwise_.size(); k++ )
		terms.push_back( pairwise_[k]->parameters() );
	int np=0;
	for( unsigned int k=0; k<pairwise_.size(); k++ )
        np += terms[k].n_rows;
    arma::vec r( np );
	for( unsigned int k=0,i=0; k<pairwise_.size(); k++ ) {
        r.subvec( i, i + terms[k].n_rows - 1 ) = terms[k];
        i += terms[k].n_rows;
	}	
	return r;
}

void DenseCRF::setLabelCompatibilityParameters( const arma::vec & v ) {
	std::vector< int > n;
	for( unsigned int k=0; k<pairwise_.size(); k++ )
        n.push_back( pairwise_[k]->parameters().n_rows );
	int np=0;
	for( unsigned int k=0; k<pairwise_.size(); k++ )
		np += n[k];
	
	for( unsigned int k=0,i=0; k<pairwise_.size(); k++ ) {
        pairwise_[k]->setParameters( v.subvec( i, i+n[k]-1 ) );
		i += n[k];
	}	
}

arma::vec DenseCRF::kernelParameters() const {
    std::vector< arma::vec > terms;
	for( unsigned int k=0; k<pairwise_.size(); k++ )
		terms.push_back( pairwise_[k]->kernelParameters() );
	int np=0;
	for( unsigned int k=0; k<pairwise_.size(); k++ )
        np += terms[k].n_rows;
    arma::vec r( np );
	for( unsigned int k=0,i=0; k<pairwise_.size(); k++ ) {
        r.subvec( i, i + terms[k].n_rows-1 ) = terms[k];
        i += terms[k].n_rows;
	}	
	return r;
}

void DenseCRF::setKernelParameters( const arma::vec & v ) {
	std::vector< int > n;
	for( unsigned int k=0; k<pairwise_.size(); k++ )
        n.push_back( pairwise_[k]->kernelParameters().n_rows );
	int np=0;
	for( unsigned int k=0; k<pairwise_.size(); k++ )
		np += n[k];
	
	for( unsigned int k=0,i=0; k<pairwise_.size(); k++ ) {
        pairwise_[k]->setKernelParameters( v.subvec( i, i+n[k]-1 ) );
		i += n[k];
	}	
}

arma::uword DenseCRF2D::getColor( const unsigned char * c ){
    return (arma::uword)qRgb(c[0],c[1],c[2]);
}

arma::uvec DenseCRF2D::getLabelingImg( const unsigned char * im, int N, int M, ColorLabelMap& map){
    arma::uvec res(N);
    //printf("%d %d %d \n",im[0],im[1],im[2]);
    ColorLabelMap tmp;
    for( int k=0; k<N; k++ ){
        // Map the color to a label
        arma::uword c = getColor( im + 3*k );
        arma::uword i;
        if(tmp.end()==tmp.find(c))
        {
            i = tmp.size()+1;
            tmp.insert(c,i);
            map.insert(i,c);
        }else{
            i = tmp[c];
        }
        res[k] = c?i:-1;
    }
    return res;
}
