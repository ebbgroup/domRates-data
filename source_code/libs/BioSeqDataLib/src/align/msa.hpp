#ifndef SRC_ALIGN_MSA_HPP_
#define SRC_ALIGN_MSA_HPP_


#include "nw_gotoh.hpp"
#include "../phylogeny/PhylogeneticTree.hpp"
#include "../DomainModule.hpp"

namespace BioSeqDataLib
{


template<typename DomainSet, typename SimMatrix, typename DataType, typename F>
void
initGotohMatrix(const DomainSet &set1, const DomainSet &set2, SimMatrix &simMat, F isGap, MatrixStack<3, std::pair<DataType, char> > &matrices)
{
	size_t len1=set1.begin().size();
	size_t len2=set2.begin().size();
	size_t itEnd1 = set1.end();
	size_t itEnd2 = set2.end();
	size_t nCombination = set1.size() * set2.size();
	matrices.resize(len1+1, len2+1);
	auto &matrix = matrices[0];
	for (size_t i=0; i<len1; ++i)
	{
		for (size_t j=0; j<len2; ++j)
		{
			DataType score = 0;
			for (auto it1=set1.begin(); it1!=itEnd1; ++it1)
			{
				for (auto it2=set2.begin(); it2!=itEnd2; ++it2)
				{
					if ((!isGap(it1->second[i])) && (!isGap(it2->second[j])))
						score += simMat.val(it1->second[i],it2->second[j]);
				}
			}
			matrix[i+1][j+1].first = score/nCombination;
		}
	}
}

/**
 * \brief Inserts gap domains into the darrangement.
 * @param domSet     The domain arrangement set
 * @param editString The edit string
 */
template <typename D>
void
insertGaps(DomainArrangementSet<D> &domSet, const std::string &editString)
{
	D emptyDomain;
	for (auto &pair : domSet)
	{
		auto &arrangement = pair.second;
		auto it = arrangement.end();
		for (size_t i = 0; i< editString.size(); ++i)
		{
			if (editString[i] == '-')
				it = arrangement.emplace(it, emptyDomain);
			else
				--it;
		}
	}
}

template <typename D, typename F, typename M, typename G>
void
progressive_align(PhylogeneticTree<std::pair<size_t, size_t> > &tree, D &dataSet, const M &simMat, F isGap, G gop, G gep)
{
	MatrixStack<3, std::pair<float, char> > matrices;
	std::vector<D> full;
	full.resize(dataSet.size());
	size_t pos = 0;
	auto itEnd = tree.postorderEnd();
	for (auto it = tree.postorderBegin(); it != itEnd; ++it)
	{
		if (it->isLeaf())
		{
			full[pos].emplace_back(dataSet[it->name]);
			if (it->childNum() == 0)
				it->parent()->data.first = pos++;
			else
				it->parent()->data.second = pos++;
		}
		else
		{
			auto &daSet1 = full[it->data.first];
			auto &daSet2 = full[it->data.second];
			auto dim1 = daSet1.begin()->size();
			auto dim2 = daSet2.begin()->size();
			std::string editString1, editString2;
			initGotohMatrix(daSet1, daSet2, simMat, isGap, matrices);
			runGotoh(matrices, dim1, dim2, gop, gep);
			gotohTraceback(matrices, dim1, dim2, editString1, editString2);
			insertGaps(daSet1, editString1);
			insertGaps(daSet2, editString2);

			for (auto &pair : daSet1)
				daSet1.emplace(std::move(pair.fist), std::move(pair.second));

			if (!it->isRoot())
			{
				if (it->childNum() == 0)
					it->parent()->data.first = it->data.first;
				else
					it->parent()->data.second = it->data.first;
			}
		}
	}
}

}

#endif
