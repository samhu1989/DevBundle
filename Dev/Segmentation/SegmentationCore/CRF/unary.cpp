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
#include "unary.h"
#include <assert.h>

UnaryEnergy::~UnaryEnergy() {
}
arma::vec UnaryEnergy::parameters() const {
    return arma::vec();
}
void UnaryEnergy::setParameters( const arma::vec & v ) {
}
arma::vec UnaryEnergy::gradient( const arma::mat & b ) const {
    return arma::vec();
}
ConstUnaryEnergy::ConstUnaryEnergy( const arma::mat & u ):unary_(u) {
}
arma::mat ConstUnaryEnergy::get() const {
	return unary_;
}

LogisticUnaryEnergy::LogisticUnaryEnergy( const arma::mat & L, const arma::mat & f ):L_(L),f_(f) {
}
arma::mat LogisticUnaryEnergy::get() const {
	return L_*f_;
}
arma::vec LogisticUnaryEnergy::parameters() const {
    return arma::vec((double*)L_.memptr(),L_.size(), true,true );
}
void LogisticUnaryEnergy::setParameters( const arma::vec & v ) {
    assert( v.n_rows == L_.size() );
    L_ = arma::mat((double*)v.memptr(),L_.n_rows,L_.n_cols,true,true);
}
arma::vec LogisticUnaryEnergy::gradient( const arma::mat & b ) const {
    arma::mat g = b*f_.t();
    return arma::vec((double*)g.memptr(),g.size(),true,true);
}
