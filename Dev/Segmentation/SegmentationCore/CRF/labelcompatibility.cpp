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
#include "labelcompatibility.h"
#include <assert.h>
LabelCompatibility::~LabelCompatibility() {
}

void LabelCompatibility::applyTranspose( arma::mat & out, const arma::mat & Q ) const {
	apply( out, Q );
}

arma::vec LabelCompatibility::parameters() const {
    return arma::vec();
}

void LabelCompatibility::setParameters( const arma::vec & v ) {
}

arma::vec LabelCompatibility::gradient( const arma::mat & b, const arma::mat & Q ) const {
    return arma::vec();
}


PottsCompatibility::PottsCompatibility( float weight ): w_(weight) {
}

void PottsCompatibility::apply( arma::mat & out, const arma::mat & Q ) const {
	out = -w_*Q;
}

arma::vec PottsCompatibility::parameters() const {
    arma::vec r(1);
	r[0] = w_;
	return r;
}

void PottsCompatibility::setParameters( const arma::vec & v ) {
	w_ = v[0];
}

arma::vec PottsCompatibility::gradient( const arma::mat & b, const arma::mat & Q ) const {
    arma::vec r(1);
    r[0] = - arma::accu( b % Q );
	return r;
}


DiagonalCompatibility::DiagonalCompatibility( const arma::vec & v ): w_(v) {
}

void DiagonalCompatibility::apply( arma::mat & out, const arma::mat & Q ) const {
    assert( w_.n_rows == Q.n_rows );
    out = arma::diagmat(w_)*Q;
}

arma::vec DiagonalCompatibility::parameters() const {
	return w_;
}

void DiagonalCompatibility::setParameters( const arma::vec & v ) {
	w_ = v;
}

arma::vec DiagonalCompatibility::gradient( const arma::mat & b, const arma::mat & Q ) const {
    return arma::sum( b%Q ,1);
}

MatrixCompatibility::MatrixCompatibility( const arma::mat & m ): w_(0.5*(m + m.t())) {
    assert( m.n_cols == m.n_rows );
}

void MatrixCompatibility::apply( arma::mat & out, const arma::mat & Q ) const {
	out = w_*Q;
}

void MatrixCompatibility::applyTranspose( arma::mat & out, const arma::mat & Q ) const {
    out = w_.t()*Q;
}

arma::vec MatrixCompatibility::parameters() const {
    arma::vec r( w_.n_cols*(w_.n_rows+1)/2 );
    for( int i=0,k=0; i<w_.n_cols; i++ )
        for( int j=i; j<w_.n_rows; j++, k++ )
			r[k] = w_(i,j);
	return r;
}

void MatrixCompatibility::setParameters( const arma::vec & v ) {
    assert( v.n_rows == w_.n_cols*(w_.n_rows+1)/2 );
    for( int i=0,k=0; i<w_.n_cols; i++ )
        for( int j=i; j<w_.n_rows; j++, k++ )
			w_(j,i) = w_(i,j) = v[k];
}

arma::vec MatrixCompatibility::gradient( const arma::mat & b, const arma::mat & Q ) const {
    arma::mat g = b * Q.t();
    arma::vec r( w_.n_cols*(w_.n_rows+1)/2 );
    for( int i=0,k=0; i<g.n_cols; i++ )
        for( int j=i; j<g.n_rows; j++, k++ )
			r[k] = g(i,j) + (i!=j?g(j,i):0.f);
	return r;
}
	
