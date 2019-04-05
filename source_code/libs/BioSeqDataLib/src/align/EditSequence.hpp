#ifndef EDITSEQUENCE_HPP
#define EDITSEQUENCE_HPP

#include <string>
#include <vector>

#include "../domain/DomainArrangement.hpp"
#include "../sequence/Sequence.hpp"



namespace BioSeqDataLib
{

/**
 * \class EditSequence
 * \brief Simple class to store alignment results.
 * 
 * The edit sequence eS1 and eS2 contain the alignment positions. Gaps are denoted as one, else the number represents the index in the sequence.
 * 
 * 0 -1 1 2  (sequence ABC) results in A-BC
 */
class EditSequence
{
    public:
    std::vector<long int> eS1; //!< Stores the edit string of sequence 1
    std::vector<long int> eS2; //!< Stores the edit string of sequence 2
    size_t start1; //!<  Start of the alignment in sequence 1
    size_t end1;   //!<  End of the alignment in sequence 1
    size_t start2; //!<  Start of the alignment in sequence 2
    size_t end2;   //!<  End of the alignment in sequence 2

    /**
     * @brief Clears the EditSequence
     */
    void
    clear()
    {
        eS1.clear();
        eS2.clear();
        start1=end1=start2=end2=0;
    }

    /**
     * @brief Returns the size of the alignment.
     * 
     * @return size_t The size of the alignment.
     */
    size_t
    size() const
    {
        return eS1.size();
    }

    /**
     * @brief Returns the the actual alignment of the two given domain arrangement
     * 
     * @return std::pair<std::string, std::string> A tuple containing the two alignment strings
     */
    template<typename D>
    std::pair<std::string, std::string>
    alnStrings(DomainArrangement<D> &da1, DomainArrangement<D> &da2)
    {
        std::string alnString1, alnString2;
        for (size_t i = 0; i<eS1.size(); ++i)
        {
            if (i>0)
            {
                alnString1 += " ";
                alnString2 += " ";
            }
            if (eS1[i] == -1)
            {
                alnString1 += std::string(da2[eS2[i]].accession().size(), '*');
                alnString2 += da2[eS2[i]].accession();
            }
            else
            {
                if (eS2[i] == -1)
                {
                    alnString1 += da1[eS1[i]].accession();
                    alnString2 += std::string(da1[eS1[i]].accession().size(), '*');
                }
                else
                {
                    alnString1 += da1[eS1[i]].accession();
                    alnString2 += da2[eS2[i]].accession();
                }
            }
        }
        return make_pair(std::move(alnString1), std::move(alnString2));
    }


};

} // BSDL namespace

#endif /*  EDITSEQUENCE  */