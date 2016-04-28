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
#include "pairwise.h"
#include <iostream>

Kernel::~Kernel() {
}
class DenseKernel: public Kernel {
protected:
	NormalizationType ntype_;
	KernelType ktype_;
	Permutohedral lattice_;
    arma::vec norm_;
    arma::mat f_;
    arma::mat parameters_;
    void initLattice( const arma::mat & f ) {
        const int N = f.n_cols;
        lattice_.init( arma::conv_to<arma::fmat>::from(f) );
		
        norm_ = arma::conv_to<arma::mat>::from(lattice_.compute( arma::fvec(N,arma::fill::ones).t() ).t());
		
		if ( ntype_ == NO_NORMALIZATION ) {
			float mean_norm = 0;
			for ( int i=0; i<N; i++ )
				mean_norm += norm_[i];
			mean_norm = N / mean_norm;
			for ( int i=0; i<N; i++ )
				norm_[i] = mean_norm;
		}
		else if ( ntype_ == NORMALIZE_SYMMETRIC ) {
			for ( int i=0; i<N; i++ )
				norm_[i] = 1.0 / sqrt(norm_[i]+1e-20);
		}
		else {
			for ( int i=0; i<N; i++ )
				norm_[i] = 1.0 / (norm_[i]+1e-20);
		}
	}
    void filter( arma::mat & out, const arma::mat & in, bool transpose ) const {
		// Read in the values
		if( ntype_ == NORMALIZE_SYMMETRIC || (ntype_ == NORMALIZE_BEFORE && !transpose) || (ntype_ == NORMALIZE_AFTER && transpose))
            out = in*arma::diagmat(norm_);
		else
			out = in;
		// Filter
        arma::fmat fout = arma::conv_to<arma::fmat>::from(out);
		if( transpose )
            lattice_.compute( fout, fout, true );
		else
            lattice_.compute( fout, fout );
        // lattice_.compute( out.data(), out.data(), out.rows() );
        out = arma::conv_to<arma::mat>::from(fout);
		// Normalize again
		if( ntype_ == NORMALIZE_SYMMETRIC || (ntype_ == NORMALIZE_BEFORE && transpose) || (ntype_ == NORMALIZE_AFTER && !transpose))
        out = out*arma::diagmat(norm_);
	}
	// Compute d/df a^T*K*b
    arma::mat kernelGradient( const arma::mat & a, const arma::mat & b ) const {
        arma::fmat g = arma::fmat(f_.n_rows,f_.n_cols,arma::fill::zeros);
        arma::fmat fa = arma::conv_to<arma::fmat>::from(a);
        arma::fmat fb = arma::conv_to<arma::fmat>::from(b);
        lattice_.gradient( (float*)g.memptr(), (float*)fa.memptr(), (float*)fb.memptr(), a.n_rows );
        return arma::conv_to<arma::mat>::from(g);
	}
    arma::mat kernelGradient( const arma::fmat & a, const arma::fmat & b ) const {
        arma::fmat g = arma::fmat(f_.n_rows,f_.n_cols,arma::fill::zeros);
        lattice_.gradient( (float*)g.memptr(), (float*)a.memptr(), (float*)b.memptr(), a.n_rows );
        return arma::conv_to<arma::mat>::from(g);
    }
    arma::mat featureGradient( const arma::mat & a, const arma::mat & b ) const {
        arma::fmat fa = arma::conv_to<arma::fmat>::from(a);
        arma::fmat fb = arma::conv_to<arma::fmat>::from(b);
        arma::fvec fnorm = arma::conv_to<arma::fvec>::from(norm_);
        arma::fmat ones = arma::fmat( a.n_rows, a.n_cols, arma::fill::ones );
        arma::mat a_n = a*arma::diagmat(norm_);
        arma::mat b_n = b*arma::diagmat(norm_);
        if (ntype_ == NO_NORMALIZATION )
			return kernelGradient( a, b );
		else if (ntype_ == NORMALIZE_SYMMETRIC ) {
            arma::fmat ffa = lattice_.compute( fa*arma::diagmat(fnorm), true );
            arma::fmat ffb = lattice_.compute( fb*arma::diagmat(fnorm) );
            arma::fvec norm3 = fnorm%fnorm%fnorm;
            arma::mat r1 = kernelGradient( 0.5*( fa%ffb + ffa%fb )*arma::diagmat(norm3), ones );
            arma::mat r2 = kernelGradient( a_n, b_n );
            return r2 - r1;
		}
		else if (ntype_ == NORMALIZE_AFTER ) {
            arma::fmat ffb = lattice_.compute( fb );
            arma::fvec norm2 = fnorm%fnorm;
            arma::mat r = kernelGradient( ( fa%ffb )*arma::diagmat(norm2), ones );
            return - r + kernelGradient( a_n, b );
		}
		else /*if (ntype_ == NORMALIZE_BEFORE )*/ {
            arma::fmat ffa = lattice_.compute( fa, true );
            arma::fmat norm2 = fnorm%fnorm;
            arma::mat r = kernelGradient( (ffa%fb)*arma::diagmat(norm2), ones );
            return -r+kernelGradient( a, b_n );
		}
	}
public:
    DenseKernel(const arma::mat & f, KernelType ktype, NormalizationType ntype):f_(f), ktype_(ktype), ntype_(ntype) {
		if (ktype_ == DIAG_KERNEL)
            parameters_ = arma::vec( f.n_rows, arma::fill::ones );
		else if( ktype == FULL_KERNEL )
            parameters_ = arma::mat( f.n_rows, f.n_rows , arma::fill::eye);
		initLattice( f );
	}
    virtual void apply( arma::mat & out, const arma::mat & Q ) const {
		filter( out, Q, false );
	}
    virtual void applyTranspose( arma::mat & out, const arma::mat & Q ) const {
		filter( out, Q, true );
	}
    virtual arma::vec parameters() const {
		if (ktype_ == CONST_KERNEL)
            return arma::vec();
		else if (ktype_ == DIAG_KERNEL)
			return parameters_;
		else {
            arma::mat p = parameters_;
            p.reshape( p.n_cols*p.n_rows, 1 );
			return p;
		}
	}
    virtual void setParameters( const arma::vec & p ) {
		if (ktype_ == DIAG_KERNEL) {
			parameters_ = p;
            initLattice( arma::diagmat(p) * f_ );
		}
		else if (ktype_ == FULL_KERNEL) {
            arma::mat tmp = p;
            tmp.reshape( parameters_.n_rows, parameters_.n_cols );
			parameters_ = tmp;
			initLattice( tmp * f_ );
		}
	}
    virtual arma::vec gradient( const arma::mat & a, const arma::mat & b ) const {
		if (ktype_ == CONST_KERNEL)
            return arma::vec();
        arma::mat fg = featureGradient( a, b );
		if (ktype_ == DIAG_KERNEL)
            return arma::sum(f_%fg,1);
		else {
            arma::mat p = fg*f_.t();
            p.reshape( p.n_cols*p.n_rows, 1 );
			return p;
		}
	}
};

PairwisePotential::~PairwisePotential(){
	delete compatibility_;
	delete kernel_;
}
PairwisePotential::PairwisePotential(const arma::mat & features, LabelCompatibility * compatibility, KernelType ktype, NormalizationType ntype) : compatibility_(compatibility) {
	kernel_ = new DenseKernel( features, ktype, ntype );
}
void PairwisePotential::apply(arma::mat & out, const arma::mat & Q) const {
	kernel_->apply( out, Q );
	// Apply the compatibility
	compatibility_->apply( out, out );
}
void PairwisePotential::applyTranspose(arma::mat & out, const arma::mat & Q) const {
	kernel_->applyTranspose( out, Q );
	// Apply the compatibility
	compatibility_->applyTranspose( out, out );
}
arma::vec PairwisePotential::parameters() const {
	return compatibility_->parameters();
}
void PairwisePotential::setParameters( const arma::vec & v ) {
	compatibility_->setParameters( v );
}
arma::vec PairwisePotential::gradient( const arma::mat & b, const arma::mat & Q ) const {
    arma::mat filtered_Q = arma::mat(Q.n_rows,Q.n_cols,arma::fill::zeros);;
	// You could reuse the filtered_b from applyTranspose
	kernel_->apply( filtered_Q, Q );
	return compatibility_->gradient(b,filtered_Q);
}
arma::vec PairwisePotential::kernelParameters() const {
	return kernel_->parameters();
}
void PairwisePotential::setKernelParameters( const arma::vec & v ) {
	kernel_->setParameters( v );
}
arma::vec PairwisePotential::kernelGradient( const arma::mat & b, const arma::mat & Q ) const {
    arma::mat lbl_Q = arma::mat(Q.n_rows,Q.n_cols,arma::fill::zeros);
	// You could reuse the filtered_b from applyTranspose
	compatibility_->apply( lbl_Q, Q );
	return kernel_->gradient(b,lbl_Q);
}
