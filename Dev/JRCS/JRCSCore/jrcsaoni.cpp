#include "jrcsaoni.h"
namespace JRCS {
JRCSAONI::JRCSAONI():JRCSBase()
{

}
void JRCSAONI::computeOnce()
{
    //reset sum
    xv_sum_.fill(0.0);
    xn_sum_.fill(0.0);
    xc_sum_.fill(0.0);
    var_sum.fill(0.0);
    alpha_sum.fill(0.0);
    alpha_sumij.fill(0.0);

    //update alpha

    //update RT

    //update X

    //fix the X center position

}
}
