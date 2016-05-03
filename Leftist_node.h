/*****************************************
 * Author: Davan Basran
 *****************************************/

#ifndef LEFTIST_NODE_H
#define LEFTIST_NODE_H

#include <algorithm>

#ifndef nullptr
#define nullptr 0
#endif

/*
class Leftist_node
Defines the nodes which contain the elements to be stored in the heap 

Member fields:
	Type element; // The element stored in the node
	Leftist_node *left_tree; // a reference to the left child, or root of the left subtree
	Leftist_node *right_tree; // a reference to the right child, or root of the right subtree
	int heap_null_path_length; // returns the null-path-length of the heap

Accessors:
	Type retrieve() const;
	returns the element stored in the node

	bool empty() const;
	returns true if the current node is null

	Leftist_node *left() const;
	retuns the left subtree

	Leftist_node *right() const;
	returns the right subtree

	int null_path_length() const;
	returns the null path length of the heap

Constuctors:
	Leftist_node<Type>::Leftist_node( Type const &obj )
	Initializes the root node and sets the left and right subtrees not be null
*/
template <typename Type>
class Leftist_node {
	private:
		Type element;
		Leftist_node *left_tree;
		Leftist_node *right_tree;
		int heap_null_path_length;

	public:
		Leftist_node( Type const & );

		Type retrieve() const;
		bool empty() const;
		Leftist_node *left() const;
		Leftist_node *right() const;
		int count( Type const & ) const;
		int null_path_length() const;

		void push( Leftist_node *, Leftist_node *& );
		void clear();
};

/*
See class definition for details
*/
template <typename Type>
Leftist_node<Type>::Leftist_node( Type const &obj ):
element( obj ),
left_tree( nullptr ),
right_tree( nullptr ),
heap_null_path_length( 0 ) {
	// does nothing
}

/*
See class definition for details
*/
template <typename Type>
bool Leftist_node<Type>::empty() const {
	return ( this == nullptr );
}

/*
See class definition for details
*/
template <typename Type>
Type Leftist_node<Type>::retrieve() const {
	return element;
}

/*
See class definition for details
*/
template <typename Type>
Leftist_node<Type> *Leftist_node<Type>::left() const {
	return left_tree;
}

/*
See class definition for details
*/
template <typename Type>
Leftist_node<Type> *Leftist_node<Type>::right() const {
	return right_tree;
}

/*
See class definition for details
*/
template <typename Type>
int Leftist_node<Type>::null_path_length() const {
	if (this == nullptr) {
		return -1;
	}
	else {
		return heap_null_path_length;
	}
}

/*
int Leftist_node<Type>::count(const Type &obj) const
Counts the number of nodes in the heap contain the obj which is passed in

Param 1:
	An object of type Type which is the obj we compare elements to 

The function has no exception handling, pre or post conditions

As long as the left or right subtree is not nullptr
	Recursivley call count on them and compare the obj to the current node, incrementing the count if they are equal
	The exit case is when both subtrees are null,those skipping the condtional statements and causing the function to return the count value
*/
template <typename Type>
int Leftist_node<Type>::count(const Type &obj) const {
	int counter = 0;

	if (empty())
		return 0;

	if (retrieve() == obj) {
		counter++;
	}
	if(left() != nullptr) counter += left()->count(obj);
	if(right() != nullptr) counter += right()->count(obj);

	return counter;
}

/*
void Leftist_node<Type>::push(Leftist_node *new_Heap, Leftist_node *&ptr_to_this)
Recursive function to push a node onto the heap

Param 1:
	A pointer to a leftist_node object. The heap that we wish to push onto another heap
Param 2:
	A pointer to a leftist_node object. The heap we are pushing an object into

Special Cases:
	If the new heap we want to push is empty then do nothing
	If the current node is empty than attach the new heap at the location of the current node

Case 1:
	If the current nodes's element is less than the element of the current node of the new heap:
		Push the new heap into the right tree
		The npl is now the smaller npl betweent left and right subtrees +1
		If the npl of the right subtree is greater than the left npl then swap trees to make it leftist
Case 2:
	The current node's element of the new heap is less than the current node's element
	Make the current node, the current node of the new heap
	push the heap into the new heap

The return has is when the new heap is empty (all elements have been pushed) thus return
*/
template <typename Type>
void Leftist_node<Type>::push(Leftist_node *new_Heap, Leftist_node *&ptr_to_this) {
	if (new_Heap == nullptr) {
		return;
	}
	if (empty()) {
		ptr_to_this = new_Heap;
		return;
	}
	if (retrieve() <= new_Heap->retrieve()) {
		if (right() == nullptr) { right_tree = new_Heap; }
		else { right_tree->push(new_Heap, right_tree); }
		heap_null_path_length = std::min(left()->null_path_length(), right()->null_path_length()) +1;
		if (left()->null_path_length() < right()->null_path_length()) {
			std::swap(left_tree, right_tree);
		}
	}
	else {
		ptr_to_this = new_Heap;
		new_Heap->push(this, new_Heap);
	}
}

/*
void Leftist_node<Type>::clear()
Recursive function which deletes all nodes the the heap

The function does not have exception handling, pre or post conditions

If the current node is null then we return (this is the exiting condition for the recursion)
Otherwise recursivley call clear on both subtrees and delete the current node
*/
template <typename Type>
void Leftist_node<Type>::clear() {
	if (this == nullptr) {
		return;
	}
	else {
		left()->clear();
		right()->clear();
		delete this;
	}
}

#endif
