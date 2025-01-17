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
#include "objective.h"
#include <cassert>
#include <algorithm>
#include <cmath>

/**** Learning Objectives ****/
ObjectiveFunction::~ObjectiveFunction(){
}
LogLikelihood::LogLikelihood( const arma::uvec & gt, double robust ):gt_( gt ),robust_(robust){
}
double LogLikelihood::evaluate( arma::mat & d_mul_Q, const arma::mat & Q ) const {
    assert( gt_.n_rows == Q.n_cols );
    const int N = Q.n_cols, M = Q.n_rows;
	double r = 0;
	d_mul_Q = 0*Q;
	for( int i=0; i<N; i++ ) 
		if( 0 <= gt_[i] && gt_[i] < M ) {
            double QQ = std::max( Q(gt_[i],i) + robust_, double(1e-20f) );
			// Make it negative since it's a 
            r += std::log(QQ) / N;
			d_mul_Q(gt_[i],i) += Q(gt_[i],i) / QQ / N;
		}
	return r;
}
Hamming::Hamming( const arma::uvec & gt, double class_weight_pow ):gt_( gt ){
    int M=0,N=gt.n_rows;
	for( int i=0; i<N; i++ )
		if( gt[i] >= M )
			M = gt[i]+1;
    arma::vec cnt(M,arma::fill::zeros);
	for( int i=0; i<N; i++ )
		if( gt[i] >= 0 )
			cnt[gt[i]] += 1;
    class_weight_ = cnt / arma::accu(cnt);
    class_weight_ = arma::pow( class_weight_ , - class_weight_pow );
    class_weight_ = class_weight_ / arma::accu( cnt % class_weight_ );
}
Hamming::Hamming( const arma::uvec & gt, const arma::mat & w ):gt_( gt ),class_weight_(w){
}
double Hamming::evaluate( arma::mat & d_mul_Q, const arma::mat & Q ) const {
    assert( gt_.n_rows == Q.n_cols );
    const int N = Q.n_cols, M = Q.n_rows;
	double r = 0;
	d_mul_Q = 0*Q;
	for( int i=0; i<N; i++ ) 
		if( 0 <= gt_[i] && gt_[i] < M ) {
			float QQ = class_weight_[ gt_[i] ] * Q(gt_[i],i);
			// Make it negative since it's a 
			r += QQ;
			d_mul_Q(gt_[i],i) += QQ;
		}
	return r;
}
IntersectionOverUnion::IntersectionOverUnion( const arma::uvec & gt ):gt_( gt ){
}
double IntersectionOverUnion::evaluate( arma::mat & d_mul_Q, const arma::mat & Q ) const {
    assert( gt_.n_rows == Q.n_cols );
    const int N = Q.n_cols, M = Q.n_rows;
	d_mul_Q = 0*Q;
	
    arma::vec in(M), un(M);
	in.fill(0.f);
	un.fill(1e-20);
	for( int i=0; i<N; i++ ) {
		if( 0 <= gt_[i] && gt_[i] < M ) {
			in[ gt_[i] ] += Q(gt_[i],i);
			un[ gt_[i] ] += 1;
			for( int l=0; l<M; l++ ) 
				if( l!=gt_[i] )
					un[ l ] += Q(l,i);
		}
	}
	for( int i=0; i<N; i++ )
		if( 0 <= gt_[i] && gt_[i] < M ) {
			for( int l=0; l<M; l++ ) 
				if( l==gt_[i] )
					d_mul_Q(l,i) = Q(l,i) / (un[l]*M);
				else
					d_mul_Q(l,i) = - Q(l,i) * in[l] / ( un[l] * un[l] * M);
		}
    return arma::accu(in%un)/M;
}

