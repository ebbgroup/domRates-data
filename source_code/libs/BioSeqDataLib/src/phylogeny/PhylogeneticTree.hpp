/*
 * PhylogenyTree.hpp
 *
 *  Created on: 10 Aug 2014
 *      Author: Carsten Kemena
 *	 Copyright: 2014
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
 * \file PhylogeneticTree.hpp
 * \brief Header containing the PhylogeneticTree.
 */
#ifndef PHYLOGENYTREE_H_
#define PHYLOGENYTREE_H_

#include <stack>
#include <vector>
#include <fstream>

#include "Tree.hpp"
#include "../utility/stringHelpers.hpp"
#include "../external/Input.hpp"


namespace BioSeqDataLib
{

/** \addtogroup PhyloGroup
 *  @{
 */

/**
 * \brief Node of a phylogenetic tree.
 */
template<typename DataType>
class TreeNodePhylo : public TreeNodeBasic<TreeNodePhylo<DataType> >
{
public:
	unsigned int id;
	std::string name;
	float edgeLength;
	float bootStrap;
	DataType data;
};


/**
 * \brief A class to represent a phylogenetic tree.
 */
template<typename DataType>
class PhylogeneticTree : public Tree<TreeNodePhylo<DataType> >
{
private:
	bool isRooted_;

	std::string readNexus_(AlgorithmPack::Input &inS);
public:

	/**
	 * \brief Standard constructor
	 */
	PhylogeneticTree();

	/**
	 * \brief Standard destructor
	 */
	virtual ~PhylogeneticTree();

	/**
	 * \brief Calculates a phylogenetic tree using the neighbour-joining algorithm.
	 * @param distMat The distance matrix.
	 * @param names The names to use.
	 */
	void nj(Matrix<float> &distMat, const std::vector<std::string> &names);

	/**
	 * \brief Calculates a phylogenetic tree using the UPGMA algorithm.
	 * @param distMat The distance matrix.
	 * @param names The names to use.
	 */
	void upgma(Matrix<float> &distMat, const std::vector<std::string> &names);

	/**
	 * \brief Turns the tree into a newick string.
	 * @return The tree in string format.
	 */
	std::string str();

	/**
	 * \brief Reads a file in newick Format.
	 * @param inFile The file to read.
	 */
	void
	read(const std::string &inFile);

	/**
	 * \brief Turns a string in newick format into a tree.
	 * @param treeLine The tree line.
	 */
	void
	str2tree(const std::string &treeLine);

	/**
	 * \brief Checks if the tree is rooted.
	 * @return true if the tree is rooted, else false.
	 */
	bool
	isRooted()
	{
		return isRooted_;
	}
};


template<typename DataType>
PhylogeneticTree<DataType>::PhylogeneticTree() : isRooted_(false)
{
	// TODO Auto-generated constructor stub

}

template<typename DataType>
PhylogeneticTree<DataType>::~PhylogeneticTree()
{
	// TODO Auto-generated destructor stub
}

template<typename DataType>
void
PhylogeneticTree<DataType>::nj(Matrix<float> &distMat, const std::vector<std::string> &names)
{
	this->root_=nullptr;
	size_t nTaxa = names.size();
	size_t nToDo = nTaxa;
	std::vector<float> r(nTaxa);
	Matrix<float> mMat(nTaxa, nTaxa);
	size_t i,j, minI=0, minJ=0;
	float minVal;
	std::vector<TreeNodePhylo<DataType> *> nodes(nTaxa);
	TreeNodePhylo<DataType> *tmpNode = nullptr;
	for (i=0; i<nTaxa; ++i)
	{
		nodes[i]= new TreeNodePhylo<DataType>();
		nodes[i]->name = names[i];
		nodes[i]->id = i;
	}
	while (nToDo >2)
	{
		// calculate r
		for (i=0; i<nTaxa; ++i)
			r[i] = 0;
		for (i=0; i<nTaxa; ++i)
		{
			if (nodes[i] == nullptr)
				continue;
			for (j=i+1; j<nTaxa; ++j)
			{
				if (nodes[j] == nullptr)
					continue;
				r[i] += distMat[i][j];
				r[j] += distMat[i][j];
			}

		}
		for (i=0; i<nTaxa; ++i)
			r[i] /= (nToDo - 2);

		// calculate intermediate matrix
		minVal = FLT_MAX;
		for (i=0; i<nTaxa; ++i)
		{
			if (nodes[i] == nullptr)
				continue;
			for (j=i+1; j<nTaxa; ++j)
			{
				if (nodes[j] == nullptr)
					continue;
				mMat[i][j]=mMat[j][i] = distMat[i][j] - (r[i] + r[j]);
				if (mMat[i][j] < minVal)
				{
					minI = i;
					minJ = j;
					minVal = mMat[i][j];
				}
			}
		}
		// constuct new Node
		tmpNode = new TreeNodePhylo<DataType>();
		tmpNode->addChild(nodes[minI]);
		tmpNode->addChild(nodes[minJ]);
		minVal = distMat[minI][minJ];
		nodes[minI]->edgeLength = (minVal + r[minI] - r[minJ])/2;
		nodes[minJ]->edgeLength =  minVal - nodes[minI]->edgeLength;

		nodes[minJ] = nullptr;
		nodes[minI] = nullptr;

		// calculate new distances
		for (i=0; i<nTaxa; ++i)
		{
			if (nodes[i] == nullptr)
				continue;
			distMat[minI][i] = distMat[i][minI] = (distMat[i][minI] + distMat[i][minJ] - minVal)/2;
		}
		--nToDo;
		nodes[minI] = tmpNode;
	}
	bool first = true;
	for (i=0; i<nTaxa; ++i)
	{
		if (nodes[i] != nullptr)
		{
			first = !first;
			if (first)
				minI = i;
			else
				minJ = i;
		}
	}

	this->root_.reset(nodes[minJ]);
	//root_->parent_ = nullptr;
	nodes[minJ]->addChild(nodes[minI]);
	nodes[minI]->edgeLength=distMat[minI][minJ];
	isRooted_=false;
}

template<typename DataType>
void
PhylogeneticTree<DataType>::upgma(Matrix<float> &distMat, const std::vector<std::string> &names)
{
	this->root_=nullptr;
	size_t nTaxa = names.size();
	size_t nToDo = nTaxa;
	std::vector<float> counts(nTaxa, 1);
	size_t i,j, minI=0, minJ=0;
	float minVal;
	std::vector<TreeNodePhylo<DataType> *> nodes(nTaxa);
	TreeNodePhylo<DataType> *tmpNode = nullptr;
	for (i=0; i<nTaxa; ++i)
	{
		nodes[i]= new TreeNodePhylo<DataType>();
		nodes[i]->edgeLength = 0;
		nodes[i]->name = names[i];
		nodes[i]->id = i;
	}

	while (nToDo >2)
	{
		// calculate intermediate matrix
		minVal = FLT_MAX;
		for (i=0; i<nTaxa; ++i)
		{
			if (nodes[i] == nullptr)
				continue;
			for (j=i+1; j<nTaxa; ++j)
			{
				if (nodes[j] == nullptr)
					continue;
				if (distMat[i][j] < minVal)
				{
					minI = i;
					minJ = j;
					minVal = distMat[i][j];
				}
			}
		}
		// constuct new Node
		tmpNode = new TreeNodePhylo<DataType>();
		tmpNode->addChild(nodes[minI]);
		tmpNode->addChild(nodes[minJ]);
		tmpNode->edgeLength = minVal/2;
		nodes[minJ]->edgeLength = minVal/2-nodes[minJ]->edgeLength;
		nodes[minI]->edgeLength = minVal/2-nodes[minI]->edgeLength;
		nodes[minJ] = nullptr;
		nodes[minI] = nullptr;

		// calculate new distances
		for (i=0; i<nTaxa; ++i)
		{
			if (nodes[i] == nullptr)
				continue;
			distMat[minI][i] = distMat[i][minI] = (distMat[i][minI]*counts[minI] + distMat[i][minJ]*counts[minJ])/(counts[minI]+counts[minJ]);
		}
		--nToDo;
		counts[minI] += counts[minJ];
		nodes[minI] = tmpNode;


	}
	bool first = true;
	for (i=0; i<nTaxa; ++i)
	{
		if (nodes[i] != nullptr)
		{
			first = !first;
			if (first)
				minI = i;
			else
				minJ = i;
		}
	}
	minVal = distMat[minI][minJ];
	this->root_.reset(new TreeNodePhylo<DataType>());
	//root_->parent = nullptr;
	this->root_->addChild(nodes[minI]);
	this->root_->addChild(nodes[minJ]);
	nodes[minJ]->edgeLength = minVal/2-nodes[minJ]->edgeLength;
	nodes[minI]->edgeLength = minVal/2-nodes[minI]->edgeLength;
}

template<typename DataType>
std::string
PhylogeneticTree<DataType>::str()
{
	std::string treeString;
	std::stack<std::pair<TreeNodePhylo<DataType>*, int> > toDo;
	toDo.push(std::make_pair(this->root_.get(), 0));
	while (!toDo.empty())
	{
		TreeNodePhylo<DataType> *currentNode = toDo.top().first;
		unsigned int child_id = toDo.top().second;
		if (!currentNode->isLeaf())
		{
			if (child_id == 0)
				treeString.push_back('(');
			else if (currentNode->nChildren() == child_id)
			{
				toDo.pop();
				if (!toDo.empty())
					treeString.append(")" + currentNode->name+ ":"+std::to_string(currentNode->edgeLength));
				else
					treeString.append(")" + currentNode->name+ ";");
				continue;
			}
			else
			{
				treeString.push_back(',');
			}
			++toDo.top().second;
			toDo.push(std::make_pair(currentNode->child(child_id),0));
		}
		else
		{
			toDo.pop();
			treeString.append(currentNode->name + ":" + std::to_string(currentNode->edgeLength));
		}
	}
	return treeString;
}

template<typename DataType>
void
PhylogeneticTree<DataType>::read(const std::string &inFile)
{
	AlgorithmPack::Input inS;
	inS.open(inFile);
	std::string treeLine;
	getline(inS, treeLine);
	if (treeLine == "#NEXUS")
		treeLine = readNexus_(inS);
	inS.close();
	str2tree(treeLine);
}

template<typename DataType>
std::string
PhylogeneticTree<DataType>::readNexus_(AlgorithmPack::Input &inS)
{
	std::string line;
	while (getline(inS, line))
	{
		if ((line[0]=='b') || (line[0] == 'B'))
		{
			std::transform(line.begin(), line.end(), line.begin(), ::toupper);
			if (!line.compare(0,11,"BEGIN TREES"))
			{
				size_t pos;
				while(true)
				{
					getline(inS, line);
					pos =0;
					while (!isalpha(line[pos]))
						++pos;
					if ((toupper(line[pos])=='T') && (toupper(line[++pos])=='R') && (toupper(line[++pos])=='E') && (toupper(line[++pos])=='E'))
						break;

				}

				line = line.substr(line.find('('), std::string::npos);
				size_t len = line.size();
				pos =0;
				for (size_t i=0; i<len; ++i)
				{
					if (line[i] == '[')
						while (line[++i] != ']');
					else
						line[pos++] = line[i];
				}
				line.resize(pos);
				return line;
			}
		}
	}
	return "";
}

template<typename DataType>
void
PhylogeneticTree<DataType>::str2tree(const std::string &treeLine)
{
	TreeNodePhylo<DataType> *currentNode = new TreeNodePhylo<DataType>();
	currentNode->id = 0;
	this->root(currentNode);
	std::stack<TreeNodePhylo<DataType> *> nodeStack;
	nodeStack.push(&this->root());
	size_t pos = 0;
	char c;
	size_t id = 0;
	std::string name;
	size_t nameEnd, nameStart;
	while ((c=treeLine[++pos]) != ';')
	{
		if (c == ',')
			continue;
		// read name if exists
		double edgeLen = 1;
		if (c != '(')
		{
			nameStart = (treeLine[pos] == ')') ? pos+1 : pos;
			nameEnd = treeLine.find_first_of(",)(;:", nameStart);
			if (nameEnd != nameStart)
				name = treeLine.substr(nameStart, nameEnd-nameStart);
			else
				name = "";
			pos = nameEnd-1;

			// add edge length
			if (treeLine[nameEnd] == ':')
			{
				size_t edgeLenStart = ++nameEnd;
				size_t edgeLenEnd = treeLine.find_first_of(",)(;", nameStart);
				edgeLen = stod(treeLine.substr(edgeLenStart, edgeLenEnd-edgeLenStart));
				pos = edgeLenEnd-1;
			}

			if (c == ')')
			{
				nodeStack.top()->name = name;
				nodeStack.top()->edgeLength = edgeLen;
				nodeStack.pop();
				continue;
			}

		}

		currentNode = new TreeNodePhylo<DataType>();
		currentNode->name = name;
		currentNode->id = ++id;
		currentNode->edgeLength = edgeLen;
		nodeStack.top()->addChild(currentNode);

		if (c == '(')
			nodeStack.push(currentNode);
	}

}



/** @} */ // PhyloGroup

} /* namespace BioSeqDataLib */

#endif /* PHYLOGENYTREE_H_ */
