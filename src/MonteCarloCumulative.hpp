#ifndef MONTECARLOCUMULATIVE_HPP_
#define MONTECARLOCUMULATIVE_HPP_

#include <vector>

class MonteCarloCumulative {
public:
	MonteCarloCumulative(const std::vector<double>& cdf, int nquantiles);
	int randomSelect();
private:
	std::vector<double> m_cdf;
	std::vector<int> m_quantileIndices;
	int binarySearchCDF(double val, int lo, int hi) const;
};


#endif /* MONTECARLOCUMULATIVE_HPP_ */
