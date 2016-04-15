#include "crf2d.h"
#include "densecrf.h"

// Certainty that the groundtruth is correct
const double GT_PROB = 0.5;

// Simple classifier that is 50% certain that the annotation is correct
arma::mat computeUnary( const arma::uvec & lbl, int M ){
    const double u_energy = -log( 1.0 / M );
    const double n_energy = -log( (1.0 - GT_PROB) / (M-1) );
    const double p_energy = -log( GT_PROB );
    arma::mat r( M, lbl.n_rows );
    r.fill(u_energy);
    //printf("%d %d %d \n",im[0],im[1],im[2]);
    for( int k=0; k<lbl.n_rows; k++ ){
        // Set the energy
        if (lbl[k]>=0){
            r.col(k).fill(n_energy);
            r(lbl[k],k) = p_energy;
        }
    }
    return r;
}

CRF2D::CRF2D(QImage& img,arma::uvec& lbl):
    input_img_(img),label_(lbl)
{
    setObjectName("CRF2D");
}

void CRF2D::process(void)
{
    int M = 21 ;
    DenseCRF2D crf(input_img_.width(),input_img_.height(),M);
    arma::mat unary = computeUnary(label_,M);
    crf.setUnaryEnergy(unary);
    arma::mat diag(unary.n_rows,unary.n_rows,arma::fill::zeros);
    diag.diag().fill(-3.0);
    crf.addPairwiseGaussian( 3, 3, new MatrixCompatibility( diag ) );
    diag.diag().fill(-10.0);
    crf.addPairwiseBilateral( 80, 80, 13, 13, 13, input_img_.bits() , new MatrixCompatibility( diag ) );
    arma::mat Q = crf.startInference();
    arma::mat t1,t2;
    QString msg;
    msg = msg.sprintf("kl = %lf\n", crf.klDivergence(Q));
    emit message(msg,0);
    for( int it=0; it<5; it++ ) {
        crf.stepInference( Q, t1, t2 );
        label_ = crf.currentMap(Q);
        msg = msg.sprintf("kl = %lf\n", crf.klDivergence(Q));
        emit message(msg,0);
    }
    emit end();
}

