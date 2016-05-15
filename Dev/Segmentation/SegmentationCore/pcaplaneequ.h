#include <vector>
#include <armadillo>
class pcaplaneequ
{
public:
	pcaplaneequ(void);
	~pcaplaneequ(void);

    int push_point(arma::fvec& p);
    arma::fvec getnormal(){ return normal;}
    arma::fvec getave(){ return ave;}
	float geta(){ return a;}
	float getb(){ return b;}
	float getc(){ return c;}
	float getd(){ return d;}
	long  getn(){ return n;}
	float geteigennormal(){ return eigennormal;}
	bool existed()
	{
		if (n>=3) return true;
		else return false;
	}
	int fitting()
	{
		if (n<3) return 0;
        arma::fvec p;
        float lambda;
		for (long i=0; i<n; i++)
		{
            p = ave - cloud[i];
            lambda = arma::dot(p,normal);
            cloud[i]+=lambda*normal;
		}
		return 1;
	}
    arma::fvec getnewpoint(arma::fvec q)
	{
        float lambda;
        arma::fvec p = ave-q;
        lambda = arma::dot(p,normal);
        q += lambda*normal;
		return q;
	}

    float dist(arma::fvec q)
	{
		return fabs(getlanda(q));
	}

    double getlanda(arma::fvec q)
	{
        float lambda;
        arma::fvec p = ave - q;
        lambda = arma::dot(p,normal);
        return lambda;
	}

    arma::fvec getnewpoint(long index)
	{
		if (index<n) return cloud[index];
	}

	void clear()
	{
		cloud.clear();
        cov_.fill(0);
        ave.fill(0);
		n=0;
		a=0;
		b=0;
		c=0;
		d=0;
		return ;
	}

    float *geteval()
	{
        return (float*)eval_.memptr();
	}

	int getsize()
	{
		return cloud.size();
	}
	

private:
    std::vector<arma::fvec> cloud;
	long n;
    arma::fvec::fixed<3> ave;
    arma::fmat::fixed<3,3> cov_;
    arma::fvec::fixed<3> eval_;
    arma::fmat::fixed<3,3> evec_;
	
	float a,b,c,d;	
    arma::fvec normal;
	float eigennormal;
};

