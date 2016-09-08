#include "MonteCarloCumulative.hpp"
#include "Random.hpp"

MonteCarloCumulative::MonteCarloCumulative(const std::vector<double>& cdf, int nquantiles)
: m_cdf(cdf)
{
	double step = 1.0 / (double)nquantiles;
	for (int i = 0; i < nquantiles - 1; ++i)
		m_quantileIndices.push_back(binarySearchCDF((i + 1) * step, 0, m_cdf.size() - 1));
}

int MonteCarloCumulative::randomSelect() {
	double rnum = Random::randomDouble(0.0, 1.0);
	int index = (int)(rnum * m_quantileIndices.size()) - 1;
	int lo = 0, hi = m_cdf.size() - 1;
	if (index >= 0)
		lo = m_quantileIndices[index + 0];
	if (index + 1 < m_quantileIndices.size())
		hi = m_quantileIndices[index + 1];
	return binarySearchCDF(rnum, lo, hi);
}

int MonteCarloCumulative::binarySearchCDF(double val, int lo, int hi) const {
	while (lo != hi) {
		int mid = lo + (hi-lo)/2;
		if (val == m_cdf[mid])
			return mid;
		else if (val > m_cdf[mid])
			lo = mid+1;
		else
			hi = mid;
	}
	return lo;
}


