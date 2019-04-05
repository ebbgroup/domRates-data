/*
 * AlignmentMatrix.hpp
 *
 *  Created on: 05.02.2018
 *      Author: CarstenK
 *		 Email: c.kemena[@]uni-muenster.de
 *	 Copyright: 2018
 *
 *  This file is part of BioSeqDataLib.
 *
 *  BioSeqDataLib is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  BioSeqDataLib is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with BioSeqDataLib.  If not, see <http://www.gnu.org/licenses/>.
 */

 /**
  * \file AlignmentMatrix.hpp
  * \brief Header containing the AlignmentMatrix class.
  */
#ifndef AlignmentMatrix_hpp
#define AlignmentMatrix_hpp

#include <map>
#include <set>
#include <tuple>
#include <vector>

#include "../utility/Matrix.hpp"
#include "../utility/MatrixStack.hpp"
#include "../utility/SimilarityMatrix.hpp"
#include "../utility/LineMatrix.hpp"
#include "EditSequence.hpp"

namespace BioSeqDataLib
{

/**
 * \class AlignmentMatrix
 * \brief A class that handled the computation of pairwise alignments
 * 
 * The class handles the usage of several different alignment algorithm. Currently supported are the 
 * Needleman-Wunsch, Smith-Waterman, Gotoh and the RASPODOM algorithm.
 * 
 */
template<typename DataType, typename SimMat>
class AlignmentMatrix {
private:
    enum class Algorithm {NW, SW, Gotoh, Raspodom_NW, Unknown}; // Class to keep track of the algorithm used

    DataType gop_; // gap opening penalty
    DataType gep_; // gap extension penalty
    SimMat simMat_;  // the similarity matrix to use

    MatrixStack<3, std::pair<DataType, char> > matrices_; // matrix for the gotoh algorithm
    mutable Matrix<std::pair<DataType, char> > matrix_;   // matrix for nw/sw/raspodom algorithm

    enum Algorithm algorithm_;  // stores the algorithm used in the last run
    size_t dim1_;   // the first dimension of the matrix
    size_t dim2_;   // the second dimension of the matrix

    // raspodom / SW, stores start position of optimal alignment
    size_t best_score_x_;  // stores the ending x-coordinate of the best alignment
    size_t best_score_y_;  // stores the ending y-coordinate of the best alignment
    size_t seq1_length_;   // length of the first sequence
    size_t seq2_length_;   // length of the second sequence


    //**********************************************************
    //*                   Result Variables                     *
    //**********************************************************

    // Variables are mutable because values might be calculated only after a traceback
    mutable DataType score_; // the final score of the program
    mutable EditSequence editString_; // stores the edit string of the alignment result
    // only for raspodom
    mutable bool isCP_;            // stores if the alignment shows a circular permutation



    //**********************************************************
    //*                      Functions                         *
    //**********************************************************

    /**
     * \brief Calculates the traceback for the NW alignment
     */
    void
    traceback_nw_() const;

    /**
     * \brief Calculates the traceback for the Gotoh alignment
     */
    void
    traceback_gotoh_() const;


    /*
     * \brief Calculates the traceback for the Raspodom alignment
     */
   // void
   // traceback_raspodom_nw_() const;

    /**
     * \brief Calculates the traceback for the Smith-Waterman alignment
     */
    void
    traceback_sw_() const;

    using MatchTracker = std::set<std::pair<size_t, size_t> >;

   /* bool
    traceback_raspodom_nw_optimal_path_(MatchTracker &matches_tl, MatchTracker &matches_tr,
        MatchTracker &matches_bl, MatchTracker &matches_br) const;

    std::pair<bool, bool>
    traceback_raspodom_nw_suboptimal_path_(size_t x, size_t y, MatchTracker &matches_tl, MatchTracker &matches_tr,
        MatchTracker &matches_bl, MatchTracker &matches_br) const;
*/

    /**
     * \brief performs a traceback
     */
    void
    traceback_() const
    {
        switch (algorithm_) {
            case Algorithm::NW:
                traceback_nw_();
                break;
            case Algorithm::Gotoh:
                traceback_gotoh_();
                break;
            /*case Algorithm::Raspodom_NW:
                traceback_raspodom_nw_();
                break;*/
            case Algorithm::SW:
                traceback_sw_();
                break;
            default:
                throw std::runtime_error("unknown algorithm");
        }
    }

public:
    /**
     * Default constructor
     */
    AlignmentMatrix(): algorithm_(Algorithm::Unknown)
    {
    }

    /**
     * \brief Constructor for using homogenous gap costs (nw, sw, raspodom)
     * @param gep    Gap extension costs
     * @param simMat Match penalty costs. Either of Type SimilarityMatrix for sequences or DSM otherwise.
     */
    AlignmentMatrix(DataType gep, const SimMat &simMat) : gep_(gep), simMat_(simMat), algorithm_(Algorithm::Unknown)
    {
    }

    /**
     * \brief Constructor for using affine gap pentalties (gotoh)
     * @param gop    Gap opening costs
     * @param gep    Gep extension pentalties
     * @param simMat Similarity matrix
     */
    AlignmentMatrix(DataType gop, DataType gep, const SimMat &simMat) : gop_(gop), gep_(gep), simMat_(simMat), algorithm_(Algorithm::Unknown)
    {}

    /**
     * \brief Default destructor
     */
    ~AlignmentMatrix()
    {}

    /**
     * \brief Runs the Needleman-Wunsch algorithm on the matrix.
     * 
     * The Running time of the algorithm is O(n*m). Only the gap extension costs are used for the calculation of the alginment.
     * Set gap opening penalties are ignored.
     * @param  seq1 First sequece
     * @param  seq2 Second sequence
     */
    template<typename SeqType>
    void
    nw(const SeqType &seq1, const SeqType &seq2);

    /**
     * \brief Run the Gotoh algorithm
     * @param  seq1 First sequence
     * @param  seq2 Second sequence
     */
    template<typename SeqType>
    void
    gotoh(const SeqType &seq1, const SeqType &seq2);

    /**
     * \brief Run the Smith-Waterman algorithm
     * @param  seq1 First sequence
     * @param  seq2 Second sequence
     */
    template<typename SeqType>
    void
    sw(const SeqType &seq1, const SeqType &seq2);

    /*
     * \brief Run the RASPDOM algorithm
     * @param  seq1 First sequence
     * @param  seq2 Second sequence
     *
    template<typename SeqType>
    void
    raspodom_nw(const SeqType &seq1, const SeqType &seq2);
    */

    /**
     * \brief Sets the gap extension penalites, used as well as homogenous gap costs.
     * @param  penalty Value of the gap costs.
     */
    void
    gep(DataType penalty)
    {
        gep_ = penalty;
    }

    /**
     * @brief Returns the current gap extension penalites
     * 
     * @return The used gap extension pentalties.traceback_raspodom_
     */
    DataType
    gep() const
    {
        return gep_;
    }

    /**
     * \brief Sets gap opening penalties.
     * @param  penalty Value of the gap costs.
     */
    void
    gop(DataType penalty)
    {
        gop_ = penalty;
    }


    /**
     * @brief Returns the current gap opening penalties
     * 
     * @return The used gap openeing penalties.
     */
    DataType
    gop() const
    {
        return gop_;
    }

    /**
     * \brief Sets the similarity Matrix.
     * @param  mat The matrix to use.
     */
    void
    scoring(const SimMat &mat)
    {
        simMat_ = mat;
    }

    /**
     * \brief Returns the score of the calculated alignment.
     * @return The score of the alignment.
     */
    DataType
    score() const
    {
        //if ((editString_.eS1.empty()) && (algorithm_ == Algorithm::Raspodom_NW))
        //    traceback_();
        return score_;
    }

    EditSequence&
    result() const
    {
        if (editString_.eS1.empty())
            traceback_();
        return editString_;
    }

    
    /**
     * \brief Returns whether the RASPDOM algorithm has determined the two sequences to be a circular permutation.
     * @return true if a circular permutation was detected, otherwise false
     *
    bool
    isCP() const
    {
        if ((editString_.eS1.empty()) && (algorithm_ == Algorithm::Raspodom_NW))
            traceback_();
        return isCP_;
    }*/


};


/***************************************************
 *                    NW - Algorithm               *
 ***************************************************/

template<typename DataType, typename SimMat>
template<typename SeqType>
void
AlignmentMatrix<DataType, SimMat>::nw(const SeqType &seq1, const SeqType &seq2)
{
    isCP_ = false;
    editString_.clear();
    editString_.start1 = 0;
    editString_.end1 = seq1.size()-1;
    algorithm_ = Algorithm::NW;
    dim1_=seq1.size();
    dim2_=seq2.size();
    size_t j;
    matrix_.resize(dim1_+1, dim2_+1);
    for (size_t i=0; i<dim1_; ++i)
    {
        for (j=0; j<dim2_; ++j)
            matrix_[i+1][j+1].first = simMat_.val(seq1[i],seq2[j]);
    }

    size_t dim2 = dim2_;
    ++dim2;
    auto itPrev = matrix_[0].begin();
    auto itVert = matrix_[0].begin();
    auto it = matrix_[0].begin();
    auto itEnd = matrix_[0].begin()+dim2;
    char path;
    size_t i;
    j=0;
    DataType score;
    for (; it!=itEnd; ++it)
        it->first=(j++)*gep_;

    for (i=1; i <= dim1_; ++i)
    {
        itPrev = it = matrix_[i].begin();
        itEnd = matrix_[i].begin()+dim2;
        itVert = matrix_[i-1].begin();
        it->first = i*gep_;
        ++it;
        for (; it!=itEnd; ++it, ++itPrev)
        {
            score=it->first+itVert->first;
            path=itVert->second;
            ++itVert;

            // check gaps;
            if (itVert->first > itPrev->first)
            {
                it->first = itVert->first;
                it->second = 'i';
            }
            else
            {
                it->first = itPrev->first;
                it->second = 'j';
            }
            it->first += gep_;

            // check match
            if ((score > it->first) || ((score == it->first) && (path=='m')))
            {
                it->first = score;
                it->second = 'm';
            }
        }
    }
    score_ = matrix_[dim1_][dim2_].first;
}


template<typename DataType, typename SimMat>
void
AlignmentMatrix<DataType, SimMat>::traceback_nw_() const
{
    auto &editString1 = editString_.eS1;
    auto &editString2 = editString_.eS2;
    editString_.start1 = 0;
    editString_.end1 = matrix_.dim1()-2;
    editString_.start2 = 0;
    editString_.end2 = matrix_.dim2()-2;

    size_t dim1 = dim1_;
    size_t dim2 = dim2_;
    while ((dim1 != 0) && (dim2 != 0))
    {
        if (matrix_[dim1][dim2].second == 'm')
        {
            --dim1;
            --dim2;
            editString1.push_back(dim1);
            editString2.push_back(dim2);
        }
        else
        {
            if (matrix_[dim1][dim2].second == 'i')
            {
                --dim1;
                editString1.push_back(dim1);
                editString2.push_back(-1);

            }
            else
            {
                --dim2;
                editString1.push_back(-1);
                editString2.push_back(dim2);
            }
        }
    }

    // reached first position of a sequence, fill remaining length with gaps
    editString1.reserve(editString1.size()+dim1);
    editString2.reserve(editString2.size()+dim1);
    while (dim1 != 0)
    {
        --dim1;
        editString1.push_back(dim1);
        editString2.push_back(-1);
    }
    editString1.reserve(editString1.size()+dim2);
    editString2.reserve(editString2.size()+dim2);
    while (dim2 != 0)
    {
        --dim2;
        editString1.push_back(-1);
        editString2.push_back(dim2);
    }
    std::reverse(editString1.begin(), editString1.end());
    std::reverse(editString2.begin(), editString2.end());

}


/***************************************************
 *                   GOTOH - Algorithm             *
 ***************************************************/


template<typename DataType, typename SimMat>
template<typename SeqType>
void
AlignmentMatrix<DataType, SimMat>::gotoh(const SeqType &seq1, const SeqType &seq2)
{
    isCP_ = false;
    editString_.clear();
    algorithm_ = Algorithm::Gotoh;
    dim1_ = seq1.size();
    dim2_ = seq2.size();
    size_t dim1 = dim1_+1;
    size_t dim2 = dim2_+1;
    matrices_.resize(dim1, dim2);
    Matrix<std::pair<DataType, char> > &matM = matrices_[0];
    Matrix<std::pair<DataType, char> > &matH = matrices_[1];
    Matrix<std::pair<DataType, char> > &matV = matrices_[2];
    typedef typename std::vector<std::pair<DataType, char> >::iterator MatrixIterator;
    MatrixIterator itM=matM[0].begin();
    MatrixIterator itH=matH[0].begin();
    MatrixIterator itV=matV[0].begin();
    MatrixIterator itVertM, itPrevH, itVertV;

    DataType use_gop;
    DataType match_score;
    const DataType MINIMUM = (-1)*(std::numeric_limits<DataType>::max()-1000);
    // 0=match, 1=insert, 2=deletion

    itM->first = 0;
    itH->first = 0;
    itV->first = 0;

    for (size_t j=1; j<dim2; ++j)
    {
        (++itM)->first = (++itH)->first=(long int)j*gep_;
        (++itV)->first = MINIMUM;
    }


    for (size_t i=1; i<dim1; ++i)
    {
        itM=matM[i].begin();
        itPrevH=itH=matH[i].begin();
        itV = matV[i].begin();
        itH->first = MINIMUM;
        itM->first = itV->first=matV[i-1][0].first+gep_;
        itVertM=matM[i-1].begin();
        itVertV=matV[i-1].begin();
        ++itH;
        ++itV;
        for (size_t j=1; j<dim2; ++j)
        {
            //calculate insert value
            match_score=itVertM->first;
            ++itVertM;
            ++itVertV;
            use_gop = (j==(dim2-1))? 0 : gop_;
            if (itVertV->first > (itVertM->first +use_gop))
            {
                itV->second = 'v';
                itV->first = itVertV->first;
            }else
            {
                itV->second = 'm';
                itV->first = itVertM->first+use_gop;
            }
            itV->first += gep_;

            //calculate deletion value
            use_gop = (i==(dim1-1))? 0 : gop_;
            if (itPrevH->first > (itM->first +use_gop))
            {
                itH->second = 'h';
                itH->first = itPrevH->first;
            } else
            {
                itH->second = 'm';
                itH->first = itM->first+use_gop;
            }
            itH->first += gep_;

            //calculate match value
            ++itM;
            match_score += simMat_.val(seq1[i-1], seq2[j-1]); //itM->first;
            if (itV->first > itH->first)
            {
                itM->second = 'v';
                itM->first = itV->first;
            }
            else
            {
                itM->second = 'h';
                itM->first = itH->first;
            }

            if (match_score >= itM->first)
            {
                itM->second = 'm';
                itM->first = match_score;
            }
            ++itH;
            ++itPrevH;
            ++itV;
        }
    }
    score_ = matrices_[0][dim1_][dim2_].first;
}



template<typename DataType, typename SimMat>
void
AlignmentMatrix<DataType, SimMat>::traceback_gotoh_() const
{
    isCP_ = false;
    size_t i = dim1_;
    size_t j = dim2_;
    editString_.start1 = 0;
    editString_.end1 = matrices_.dim1()-2;
    editString_.start2 = 0;
    editString_.end2 = matrices_.dim2()-2;
    auto &editString1 = editString_.eS1;
    auto &editString2 = editString_.eS2;
    int mat = 0;
    while ((i!=0) && (j!=0))
    {
        char state = matrices_[mat][i][j].second;
        if (mat==0)
        {
            if (state=='m')
            {
                --i;
                --j;
                editString1.push_back(i);
                editString2.push_back(j);
            }
            else
            {
                if (state=='v')
                    mat=2;
                else
                    mat=1;
            }
        }
        else
        {
            if (mat==2)
            {
                --i;
                editString1.push_back(i);
                editString2.push_back(-1);
            }
            else
            {
                --j;
                editString1.push_back(-1);
                editString2.push_back(j);
            }
            if (state=='m')
                mat = 0;
        }
    }

    while (j>0)
    {
        --j;
        editString1.push_back(-1);
        editString2.push_back(j);
    }
    while (i>0)
    {
        --i;
        editString1.push_back(i);
        editString2.push_back(-1);
    }

    std::reverse(editString1.begin(), editString1.end());
    std::reverse(editString2.begin(), editString2.end());
}


/***************************************************
 *               RASPODOM - Algorithm               *
 ***************************************************/
/*
template<typename DataType, typename SimMat>
template<typename SeqType>
void
AlignmentMatrix<DataType, SimMat>::raspodom_nw(const SeqType &seq1, const SeqType &seq2)
{
    editString_.clear();
    algorithm_ = Algorithm::Raspodom_NW;
    dim1_= 2 * seq1.size();
    dim2_= 2 * seq2.size();
    seq1_length_ = seq1.size();
    seq2_length_ = seq2.size();
    size_t j;
    matrix_.resize(dim1_+1, dim2_+1);
    for (size_t i=0; i<dim1_; ++i)
    {
        for (j=0; j<dim2_; ++j)
            matrix_[i+1][j+1].first = simMat_.val(seq1[i % seq1.size()],seq2[j % seq2.size()]);
    }

    size_t dim2 = dim2_;
    ++dim2;
    auto itPrev = matrix_[0].begin();
    auto itVert = matrix_[0].begin();
    auto it = matrix_[0].begin();
    auto itEnd = matrix_[0].begin()+dim2;
    char path;
    size_t i;
    j=0;
    DataType score;
    for (; it!=itEnd; ++it)
        it->first=0;

    for (i=1; i <= dim1_; ++i)
    {
        itPrev = it = matrix_[i].begin();
        itEnd = matrix_[i].begin()+dim2;
        itVert = matrix_[i-1].begin();
        it->first = 0;
        ++it;
        for (; it!=itEnd; ++it, ++itPrev)
        {
            score=it->first+itVert->first;
            path=itVert->second;
            ++itVert;
            // check gaps;
            if (itVert->first > itPrev->first)
            {
                it->first = itVert->first;
                it->second = 'i';
            }
            else
            {
                it->first = itPrev->first;
                it->second = 'j';
            }
            it->first += gep_;
            // check match
            if ((score > it->first) || ((score == it->first) && (path=='m')))
            {
                it->first = score;
                it->second = 'm';
            }
        }
    }

    // Calculate best score in bottom left quadrant and store.
    score_ = 0;
    for (size_t i=seq1.size(); i<=dim1_; ++i)
    {
        for (size_t j=seq2.size(); j<=dim2_; ++j)
        {
            if (matrix_[i][j].first > score_)
            {
                score_ = matrix_[i][j].first;
                best_score_x_ = i;
                best_score_y_ = j;
            }
        }
    }
}


template<typename DataType, typename SimMat>
bool
AlignmentMatrix<DataType, SimMat>::traceback_raspodom_nw_optimal_path_(MatchTracker &matches_tl, MatchTracker &matches_tr,
        MatchTracker &matches_bl, MatchTracker &matches_br)  const
{
    editString_.start1 = 0;
    editString_.end1 = matrix_.dim1()-2;
    editString_.start2 = 0;
    editString_.end2 = matrix_.dim2()-2;
    auto &editString1 = editString_.eS1;
    auto &editString2 = editString_.eS2;
    size_t dim1 = best_score_x_;
    size_t dim2 = best_score_y_;
    size_t stop1 = best_score_x_ % seq1_length_;
    if (stop1 == 0)
        stop1 = seq1_length_;
    size_t stop2 = best_score_y_ % seq2_length_;
    if (stop2 == 0)
        stop2 = seq2_length_;
    while ((dim1 > stop1) || (dim2 > stop2))
    {
        std::cout << dim1 << " " << dim2 << std::endl;
        if (matrix_[dim1][dim2].second == 'm')
        {
            if ((dim1 > seq1_length_) && (dim2 > seq2_length_))
                matches_br.emplace(dim1, dim2);
            else
            {
                if (dim1 > seq1_length_)
                    matches_tr.emplace(dim1, dim2);
                else
                {
                    if (dim2 > seq2_length_)
                        matches_bl.emplace(dim1, dim2);
                    else
                        matches_tl.emplace(dim1, dim2);
                }
            }
            matrix_[dim1][dim2].second = 'd';
            --dim1;
            --dim2;
            editString1.push_back(dim1 % seq1_length_);
            editString2.push_back(dim2 % seq2_length_);
        }
        else
        {
            if (matrix_[dim1][dim2].second == 'i')
            {
                matrix_[dim1][dim2].second = 'd';
                --dim1;
                editString1.push_back(dim1 % seq1_length_);
                editString2.push_back(-1);
            }
            else
            {
                matrix_[dim1][dim2].second = 'd';
                --dim2;
                editString1.push_back(-1);
                editString2.push_back(dim2 % seq2_length_);
            }
        }
    }
    // reached first position of a sequence, fill remaining length with gaps
    editString1.reserve(editString1.size() + dim1);
    editString2.reserve(editString2.size() + dim1);
    while (dim1 > stop1)
    {
        matrix_[dim1][dim2].second = 'd';
        --dim1;
        editString1.push_back(dim1 % seq1_length_);
        editString2.push_back(-1);
    }
    editString1.reserve(editString1.size() + dim2);
    editString2.reserve(editString2.size() + dim2);
    while (dim2 > stop2)
    {
        matrix_[dim1][dim2].second = 'd';
        --dim2;
        editString1.push_back(-1);
        editString2.push_back(dim2 % seq2_length_);
    }
    std::reverse(editString1.begin(), editString1.end());
    std::reverse(editString2.begin(), editString2.end());


    if (matrix_[dim1][dim2].second == 'm')
    {
        matches_tl.emplace(dim1, dim2);
        return true;
    }
    return false;
}


template<typename DataType, typename SimMat>
std::pair<bool, bool>
AlignmentMatrix<DataType, SimMat>::traceback_raspodom_nw_suboptimal_path_(size_t x, size_t y, MatchTracker &matches_tl, MatchTracker &matches_tr,
        MatchTracker &matches_bl, MatchTracker &matches_br)  const
{
    matches_tl.clear();
    matches_tr.clear();
    matches_bl.clear();
    matches_br.clear();
    size_t dim1 = x;
    size_t dim2 = y;
    size_t stop1 = x % seq1_length_;
    if (stop1 == 0)
        stop1 = seq1_length_;
    size_t stop2 = y % seq2_length_;
    if (stop2 == 0)
        stop2 = seq2_length_;

    while ((dim1 > stop1) && (dim2 >stop2))
    {
        if (matrix_[dim1][dim2].second == 'm')
        {
            if ((dim1 > seq1_length_) && (dim2 > seq2_length_))
                matches_br.emplace(dim1, dim2);
            else
            {
                if (dim1 > seq1_length_)
                    matches_tr.emplace(dim1, dim2);
                else
                {
                    if (dim2 > seq2_length_)
                        matches_bl.emplace(dim1, dim2);
                    else
                        matches_tl.emplace(dim1, dim2);
                }
            }
            matrix_[dim1][dim2].second = 'd';
            --dim1;
            --dim2;
        }
        else
        {
            if (matrix_[dim1][dim2].second == 'i')
            {
                matrix_[dim1][dim2].second = 'd';
                --dim1;
            }
            else
            {
                if (matrix_[dim1][dim2].second == 'j')
                {
                    matrix_[dim1][dim2].second = 'd';
                    --dim2;
                }
                else
                    return std::make_pair(false, false);
            }
        }
    }
    while (dim1 > stop1)
    {
        matrix_[dim1][dim2].second = 'd';
        --dim1;
    }
    while (dim2 > stop2)
    {
        matrix_[dim1][dim2].second = 'd';
        --dim2;
    }
    if (matrix_[dim1][dim2].second == 'm')// && (matches_tl.find(std::make_pair(dim1, dim2))!= matches_tl.end()))
    {
        matches_tl.emplace(dim1, dim2);
        return std::make_pair(true, true);
    }
    return  std::make_pair(false, true);
}


template<typename DataType, typename SimMat>
void
AlignmentMatrix<DataType, SimMat>::traceback_raspodom_nw_() const
{
    MatchTracker matches_tl, matches_tr, matches_bl, matches_br;
    isCP_ = false;
    short path = 0;
    // check optimal path
    bool potential_CP = true;
    bool match = traceback_raspodom_nw_optimal_path_(matches_tl, matches_tr, matches_bl, matches_br);
        // check if three quadrants are touched
    std::cout << "OPTIMAL: "<< matches_br.size() << " " << matches_bl.size() << " " << matches_tl.size() << " " << matches_tr.size() <<"\n";
    if ((matches_tl.size() != 0) && (matches_br.size() != 0) && ((matches_bl.size() != 0) || (matches_tr.size() != 0)))
    {
        if (!match)
            potential_CP = false;
        path = (matches_bl.size() != 0) ? 1 : 2;
    }
    else
        potential_CP = false;

    std::cout << potential_CP << " p\n";

    int path_counter = 0;
    if (potential_CP)
    {
        std::set<std::tuple<DataType, size_t, size_t>, std::greater<std::tuple<DataType, size_t, size_t> > > scores;
        for (size_t i=seq1_length_; i<=dim1_; ++i)
        {
            for (size_t j=seq2_length_; j<=dim2_; ++j)
            {
                if (matrix_[i][j].second == 'm')
                    scores.emplace(matrix_[i][j].first, i, j);
            }
        }

        for (const auto &elem : matches_br)
            scores.erase(std::make_tuple(matrix_[elem.first][elem.second].first, elem.first, elem.second));

        for (const auto &elem : scores)
        {
            potential_CP = false;
            auto success = traceback_raspodom_nw_suboptimal_path_(std::get<1>(elem), std::get<2>(elem), matches_tl,  matches_tr, matches_bl, matches_br);
            std::cout << success.first << " " << success.second << std::endl;
            if (success.second)
            {
                for (const auto &elem : matches_br)
                    scores.erase(std::make_tuple(matrix_[elem.first][elem.second].first, elem.first, elem.second));
                if ((matches_tl.size() != 0) && (matches_br.size() != 0) && ((matches_bl.size() != 0) || (matches_tr.size() != 0)))
                {
                    if (success.first)
                    {
                        if (((matches_bl.size() != 0)  && (path ==2)) || ((matches_tr.size() != 0)  && (path == 1)))
                            potential_CP = true;
                        ++path_counter;
                    }
                }
                if (potential_CP)
                    isCP_ = true;
            }
        }
    }
    std::cout << path_counter << "PC\n";
    if (isCP_ && (path_counter == 1))
        isCP_ = true;
}
*/

/***************************************************
 *                    SW - Algorithm               *
 ***************************************************/

template<typename DataType, typename SimMat>
template<typename SeqType>
void
AlignmentMatrix<DataType, SimMat>::sw(const SeqType &seq1, const SeqType &seq2)
{
    isCP_ = false;
    editString_.clear();
    algorithm_ = Algorithm::SW;
    dim1_ = seq1.size();
    dim2_ = seq2.size();
    size_t j;
    matrix_.resize(dim1_+1, dim2_+1);
    for (size_t i=0; i<dim1_; ++i)
    {
        for (j=0; j<dim2_; ++j)
            matrix_[i+1][j+1].first = simMat_.val(seq1[i],seq2[j]);
    }
    size_t dim2 = dim2_ + 1;
    auto itPrev = matrix_[0].begin();
    auto itVert = matrix_[0].begin();
    auto it = matrix_[0].begin();
    auto itEnd = matrix_[0].begin()+dim2;
    char path;
    size_t i;
    DataType score;
    for (; it!=itEnd; ++it)
    {
        it->first=0;
        it->second='o';
    }

    for (i=1; i <= dim1_; ++i)
    {
        itPrev = it = matrix_[i].begin();
        itEnd = matrix_[i].begin()+dim2;
        itVert = matrix_[i-1].begin();
        it->first = 0;
        it->second = 'o';
        ++it;
        for (; it!=itEnd; ++it, ++itPrev)
        {
            score=it->first+itVert->first;
            path=itVert->second;
            ++itVert;

            // check gaps;
            if (itVert->first > itPrev->first)
            {
                it->first = itVert->first;
                it->second = 'i';
            }
            else
            {
                it->first = itPrev->first;
                it->second = 'j';
            }
            it->first += gep_;

            // check match
            if ((score > it->first) || ((score == it->first) && (path=='m')))
            {
                it->first = score;
                it->second = 'm';
            }

            if (it->first <= 0)
            {
                it->first = 0;
                it->second = 'o';
            }
        }
    }

    // set best score
    score_ = 0;
    for (size_t i = 1; i <= dim1_; ++i)
    {
        for (size_t j = 1; j <= dim2_; ++j)
        {
            if (matrix_[i][j].first > score_)
            {
                score_ = matrix_[i][j].first;
                best_score_x_ = i;
                best_score_y_ = j;
            }
        }
    }
}

template<typename DataType, typename SimMat>
void
AlignmentMatrix<DataType, SimMat>::traceback_sw_() const
{
    auto &editString1 = editString_.eS1;
    auto &editString2 = editString_.eS2;
    size_t dim1 = best_score_x_;
    size_t dim2 = best_score_y_;
    editString_.end1 = dim1-1;
    editString_.end2 = dim2-1;

    while (matrix_[dim1][dim2].second != 'o')
    {
        if (matrix_[dim1][dim2].second == 'm')
        {
            --dim1;
            --dim2;
            editString1.push_back(dim1);
            editString2.push_back(dim2);
        }
        else
        {
            if (matrix_[dim1][dim2].second == 'i')
            {
                --dim1;
                editString1.push_back(dim1);
                editString2.push_back(-1);
            }
            else
            {
                --dim2;
                editString1.push_back(-1);
                editString2.push_back(dim2);
            }
        }
    }
    editString_.start1=dim1;
    editString_.start2=dim2;
    std::reverse(editString1.begin(), editString1.end());
    std::reverse(editString2.begin(), editString2.end());
}

}

#endif /* AlignmentMatrix_hpp */
