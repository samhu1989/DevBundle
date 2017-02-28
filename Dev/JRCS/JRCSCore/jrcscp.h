#ifndef JRCSCP_H
#define JRCSCP_H
/*
 * JRCS Primitive with more constraints
 */
#include "jrcsprimitive.h"
namespace JRCS {
class JRCSCORESHARED_EXPORT CP:public Plate
{
public:
    typedef enum{
        H, // horizontal
        V  // vertical
    }T;
};
class JRCSCORESHARED_EXPORT JRCSCP:public JRCSPrimitive
{
public:
    JRCSCP();
protected:
private:
};
}
#endif // JRCSCP_H
