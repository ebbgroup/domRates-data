/*
 * Tree.hpp
 *
 *  Created on: 24 Mar 2014
 *      Author: Carsten Kemena
 *		 Email: c.kemena[@]uni-muenster.de
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
 * \file Tree.hpp
 * \brief Header containing the Tree.
 */

#ifndef TREE_HPP_
#define TREE_HPP_

#include <cfloat>
#include <iostream>
#include <fstream>
#include <memory>
#include <stack>
#include <string>
#include <stack>
#include <vector>

#include "../utility/Matrix.hpp"

namespace BioSeqDataLib
{

/** \addtogroup PhyloGroup
 *  @{
 */

/**
 * \brief Basic Tree node.
 */
template<typename DataType>
class TreeNodeBasic
{
private:
	size_t childNum_;
	DataType* parent_;
	std::vector<std::shared_ptr<DataType> > children_;

public:

	TreeNodeBasic() : childNum_(-1), parent_(nullptr)
	{}

	/**
	 * \brief Returns pointer to the child.
	 * @param i The index of the child.
	 * @return Pointer to the child.
	 */
	DataType *
	child(int i)
	{
		return children_[i].get();
	}

	/**
	 * \brief Returns pointer to the child.
	 * @param i The index of the child.
	 * @return Pointer to the child.
	 */
	const DataType *
	child(int i) const
	{
		return children_[i].get();
	}

	/**
	 * \brief Return pointer to the parent.
	 * @return Pointer to the parent.
	 */
	DataType *
	parent()
	{
		return parent_;
	}

	const DataType *
	parent() const
	{
		return parent_;
	}

	/**
	 * \brief Returns the child number.
	 * @return The child number.
	 */
	size_t
	childNum() const
	{
		return childNum_;
	}

	/**
	 * \brief Adds a new node.
	 * @param node The node to add.
	 */
	void addChild(DataType *node)
	{
		node->childNum_ = children_.size();
		node->parent_ = static_cast<DataType *>(this);
		children_.emplace_back(node);
	}


	/**
	 * \brief Returns whether the node is a leaf.
	 * @return true if the node is a leaf, else false.
	 */
	bool isLeaf() const
	{
		return children_.empty();
	}

	/**
	 * \brief Returns whether the node is root.
	 * @return true if the node is the root, else false.
	 */
	bool isRoot()
	{
		return parent_ == nullptr;
	}

	/**
	 * \brief Returns the number of children.
	 * @return The number of children.
	 */
	size_t nChildren() const
	{
		return children_.size();
	}
};


/**
 * \brief Tree class
 * \tparam The internal NodeType to use
 */
template<typename TreeNodeType>
class Tree
{
protected:
	std::shared_ptr<TreeNodeType> root_;

public:

	typedef TreeNodeType value_type;

	/**
	 * \brief Standard constructor
	 */
	Tree();

	/**
	 * \brief Standard constructor
	 */
	virtual ~Tree();

	/**
	 * \brief Returns the root
	 * @return The root
	 */
	TreeNodeType&
	root()
	{
		return *(root_.get());
	}

	const TreeNodeType&
	root() const
	{
		return *(root_.get());
	}

	/**
	 * \brief Sets a new root.
	 * @param newRoot The new root.
	 */
	void
	root(TreeNodeType *newRoot)
	{
		root_.reset(newRoot);
	}



	template<bool is_const_iterator = true>
	class preorder_iterator_helper: public std::iterator<std::bidirectional_iterator_tag, TreeNodeType>
	{
			friend class Tree<TreeNodeType>;
	private:
		typedef typename std::conditional<is_const_iterator, const TreeNodeType*, TreeNodeType*>::type IteratorDataType;
		typedef typename std::conditional<is_const_iterator, const value_type, value_type>::type IteratorValueType;
		IteratorDataType _current;
		IteratorDataType _last;
		/**
		 * \brief Constructor pointing to a given element
		 * @param p Pointer to the element.
		 */
		explicit preorder_iterator_helper(IteratorDataType p) throw():_current(p)
		{}



	public:

		/**
		 * \brief Empty constructor
		 */
		preorder_iterator_helper():_current(nullptr)
		{}

		/**
		 * \brief Copy constructor
		 * @param other The iterator that should be duplicated.
		 */
		preorder_iterator_helper(const preorder_iterator_helper<false>& other):_current(other._current)
		{}

		/**
		 *  \brief Standard destructor
		 */
		virtual ~preorder_iterator_helper()
		{}



		/**
		 * \brief Compares two iterators for equality.
		 * @param iter The other iterator.
		 * @return true if both iterators point to the same element else false.
		 */
		bool operator==(const preorder_iterator_helper & iter) const throw()
		{
			return _current == iter._current;
		}

		/**
		 * \brief Compares two iterators for inequality.
		 * @param iter The other iterator.
		 * @return true if iterators point to different elements else false.
		 */
		bool operator != (const preorder_iterator_helper & iter) const throw()
		{
			return !(*this == iter);
		}

		/**
		 * \brief Prefix increment operator.
		 * @return Current iterator.
		 */
		preorder_iterator_helper & operator++()
		{
			_last = _current;
			if (_current->isLeaf())
			{
				while (_current->parent() != nullptr)
				{
					if (_current->childNum() == (_current->parent()->nChildren()-1))
						_current = _current->parent();
					else
					{
						_current = _current->parent()->child(_current->childNum()+1);
						return *this;
					}
				}
				_current=nullptr;
			}
			else
				_current = _current->child(0);
			return *this;
		}

		/**
		 * \brief Postfix increment operator.
		 * @return Current iterator.
		 */
		preorder_iterator_helper operator++(int)
		{
			const preorder_iterator temp(*this);
			operator++();
			return temp;
		}

		/**
		 * \brief Prefix decrement operator.
		 * @return Current iterator.
		 */
		preorder_iterator_helper & operator--()
		{
			if (_current==nullptr)
				_current = _last;
			else
			{
				if (_current->childNum() == 0)
					_current = _current->parent();
				else
				{
					_current = _current->parent()->child(_current->childNum()-1);
					while (!_current->isLeaf())
						_current = _current->child(_current->nChildren()-1);
				}
			}
			return *this;
		}

		/**
		 * \brief Postfix decrement operator.
		 * @return Iterator
		 */
		preorder_iterator_helper operator--(int)
		{
			const preorder_iterator temp(*this);
			operator--();
			return temp;
		}

		/**
		 * \brief Dereferencing of iterator.
		 * @return Reference to the element pair.
		 */
		IteratorValueType & operator*()throw()
		{
			return *_current;
		}

		/**
		 * \brief Dereferencing of iterator.
		 * @return Const reference to the element pair.
		 */
		IteratorValueType& operator*() const
		{
			return *_current;
		}

		/**
		 * \brief Structure dereference operator.
		 * @return Pointer to the element pair.
		 */
		IteratorValueType* operator->() const
		{
			return _current;
		}

		friend class preorder_iterator_helper<true>;

	};

	/**
	 * \name Iterators
	 */
	/**@{*/

	/**
	 * \brief Bidirectional iterator.
	 */
	typedef preorder_iterator_helper<false> preorder_iterator;

	/**
	 * \brief const bidirectional iterator.
	 */
	typedef preorder_iterator_helper<true> const_preorder_iterator;

	/**
	 * \brief reverse bidirectional iterator.
	 */
	typedef std::reverse_iterator<preorder_iterator_helper<false> > reverse_preorder_iterator;

	/**
	 * \brief const reverse bidirectional iterator.
	 */
	typedef std::reverse_iterator<preorder_iterator_helper<true> > const_reverse_preorder_iterator;
	/**@}*/



	preorder_iterator preorderBegin() throw()
	{
		return preorder_iterator(&(this->root()));
	}

	/**
	 * \brief Returns an iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	preorder_iterator preorderEnd() throw()
	{
		preorder_iterator tmp(nullptr);
		tmp._last = &(this->root());
		while (!tmp._last->isLeaf())
		{
			tmp._last = tmp._last->child(tmp._last->nChildren()-1);
		}
		return tmp;
	}

	/**
	 * \brief Returns a const_iterator to the first element in the container.
	 * @return An iterator to the first element in the container.
	 */
	const_preorder_iterator preorderBegin() const throw()
	{
		return const_preorder_iterator(&this->root());
	}


	/**
	 * \brief Returns a const_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	const_preorder_iterator preorderEnd() const throw()
	{
		const_preorder_iterator tmp(nullptr);
		tmp._last = &(this->root());
		while (!tmp._last->isLeaf())
		{
			tmp._last = tmp._last->child(tmp._last->nChildren()-1);
		}
		return tmp;
	}



	template<bool is_const_iterator = true>
	class postorder_iterator_helper: public std::iterator<std::bidirectional_iterator_tag, TreeNodeType>
	{
		friend class Tree<TreeNodeType>;
	private:
		typedef typename std::conditional<is_const_iterator, const TreeNodeType*, TreeNodeType*>::type IteratorDataType;
		typedef typename std::conditional<is_const_iterator, const value_type, value_type>::type IteratorValueType;
		IteratorDataType _current;
		IteratorDataType _last;
		/**
		 * \brief Constructor pointing to a given element
		 * @param p Pointer to the element.
		 */
		explicit postorder_iterator_helper(IteratorDataType p) throw():_current(p)
		{}



	public:

		/**
		 * \brief Empty constructor
		 */
		postorder_iterator_helper():_current(nullptr)
		{}

		/**
		 * \brief Copy constructor
		 * @param other The iterator that should be duplicated.
		 */
		postorder_iterator_helper(const postorder_iterator_helper<false>& other):_current(other._current)
		{}

		/**
		 *  \brief Standard destructor
		 */
		virtual ~postorder_iterator_helper()
		{}



		/**
		 * \brief Compares two iterators for equality.
		 * @param iter The other iterator.
		 * @return true if both iterators point to the same element else false.
		 */
		bool operator==(const postorder_iterator_helper & iter) const throw()
		{
			return _current == iter._current;
		}

		/**
		 * \brief Compares two iterators for inequality.
		 * @param iter The other iterator.
		 * @return true if iterators point to different elements else false.
		 */
		bool operator != (const postorder_iterator_helper & iter) const throw()
		{
			return !(*this == iter);
		}

		/**
		 * \brief Prefix increment operator.
		 * @return Current iterator.
		 */
		postorder_iterator_helper & operator++()
		{

			if (_current->parent() == nullptr)
			{
				_last = _current;
				_current = nullptr;
			}
			else
			{
				if (_current->childNum() == _current->parent()->nChildren()-1)
					_current = _current->parent();
				else
				{
					_current = _current->parent()->child(_current->childNum()+1);
					while (!_current->isLeaf())
						_current = _current->child(0);
			}
			}
			return *this;
		}

		/**
		 * \brief Postfix increment operator.
		 * @return Current iterator.
		 */
		postorder_iterator_helper operator++(int)
		{
			const postorder_iterator temp(*this);
			operator++();
			return temp;
		}

		/**
		 * \brief Prefix decrement operator.
		 * @return Current iterator.
		 */
		postorder_iterator_helper & operator--()
		{
			if (_current == nullptr)
				_current = _last;
			else
			{
				if (!_current->isLeaf())
					_current = _current->child(_current->nChildren()-1);
				else
				{
					while (_current->childNum() == 0)
						_current = _current->parent();
					_current = _current->parent()->child(_current->childNum()-1);

				}
			}
			return *this;
		}

		/**
		 * \brief Postfix decrement operator.
		 * @return Iterator
		 */
		postorder_iterator_helper operator--(int)
		{
			const postorder_iterator temp(*this);
			operator--();
			return temp;
		}

		/**
		 * \brief Dereferencing of iterator.
		 * @return Reference to the element pair.
		 */
		IteratorValueType & operator*()throw()
		{
			return *_current;
		}

		/**
		 * \brief Dereferencing of iterator.
		 * @return Const reference to the element pair.
		 */
		IteratorValueType& operator*() const
		{
			return *_current;
		}

		/**
		 * \brief Structure dereference operator.
		 * @return Pointer to the element pair.
		 */
		IteratorValueType* operator->() const
		{
			return _current;
		}

		friend class postorder_iterator_helper<true>;

	};

	/**
	 * \name Iterators
	 */
	/**@{*/

	/**
	 * \brief Bidirectional iterator.
	 */
	typedef postorder_iterator_helper<false> postorder_iterator;

	/**
	 * \brief const bidirectional iterator.
	 */
	typedef postorder_iterator_helper<true> const_postorder_iterator;

	/**
	 * \brief reverse bidirectional iterator.
	 */
	typedef std::reverse_iterator<postorder_iterator_helper<false> > reverse_postorder_iterator;

	/**
	 * \brief const reverse bidirectional iterator.
	 */
	typedef std::reverse_iterator<postorder_iterator_helper<true> > const_reverse_postorder_iterator;
	/**@}*/



	postorder_iterator postorderBegin() throw()
	{
		postorder_iterator tmp(&this->root());
		while (!tmp._current->isLeaf())
			tmp._current = tmp._current->child(0);
		return tmp;
	}

	/**
	 * \brief Returns an iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	postorder_iterator postorderEnd() throw()
	{
		postorder_iterator tmp(nullptr);
		tmp._last = &(this->root());
		return tmp;
	}

	/**
	 * \brief Returns a const_iterator to the first element in the container.
	 * @return An iterator to the first element in the container.
	 */
	const_postorder_iterator postorderBegin() const throw()
	{
		const_postorder_iterator tmp(&this->root());
		tmp._current = &(this->root());
		while (!tmp._current->isLeaf())
			tmp._current = tmp._current->child(0);
		return tmp;
	}


	/**
	 * \brief Returns a const_iterator to the past-the-end element in the container.
	 * @return An iterator to the past-the-end element.
	 */
	const_postorder_iterator postorderEnd() const throw()
	{
		const_postorder_iterator tmp(nullptr);
		tmp._last = &(this->root());
		return tmp;
	}

};

template<typename TreeNodeType>
Tree<TreeNodeType>::Tree()
{
	// TODO Auto-generated destructor stub
}


template<typename TreeNodeType>
Tree<TreeNodeType>::~Tree()
{
	// TODO Auto-generated destructor stub
}

/**
 * \brief Checks if a tree is strictly bifurcating.
 * \details For a strictly bifurcating tree every inner node can just have exactly two child nodes;
 * based on the throwError parameter the function either returns boolean (true = tree is strictly bifurcating)
 * or throws at first non-bifurcating node a runtime error and quits.
 *
 * \param tree tree to be checked for bifurcation
 * \param throwError if set to true a runtime error is thrown at first non-bifurcating node instead of a false return value
 *
 * \relates Tree
 */
template <typename DataType>
bool isBifurcatingTree(const Tree<DataType> &tree, bool throwError = false) {
	bool bifur = true;

	for (auto cNode = tree.postorderBegin(); cNode != tree.postorderEnd(); ++cNode) {
		if (!cNode->isLeaf()) {
			unsigned int cnum = cNode->nChildren();
			if (cnum == 1) {
				if (throwError) {
					throw std::runtime_error("Error: Not a bifurcating tree. Node with just one child detected - Node-ID " + std::to_string(cNode->id));
				}
				else {
					bifur = false;
				}
			}
			if (cnum > 2) {
				if (throwError) {
					throw std::runtime_error("Error: Not a bifurcating tree. Node with more than two childs detected - Node-ID " + std::to_string(cNode->id));
				}
				else {
					bifur = false;
				}
			}
		}
	}

	return bifur;
}

/** @} */ // PhyloGroup

} /* namespace BioSeqDataLib */

#endif /* TREE_HPP_ */
