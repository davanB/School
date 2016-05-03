/*****************************************
 * Author: Davan Basran
 *****************************************/

#ifndef LEFTIST_HEAP_H
#define LEFTIST_HEAP_H

#ifndef nullptr
#define nullptr 0
#endif

#include "Leftist_node.h"
#include "Dynamic_queue.h"
#include <stack> // for the cout override 

/*
class Leftist_heap
Holds a heap object which contains pointers of heap_node objects

Member FIelds:
	Leftist_node *root_node // a reference to the root node of the heap
	int heap_size // contins the count of the number of nodes in the heap

Accessors:
	bool empty() const;
	return true if the heap is currently empty

	int size() const;
	returns the number of nodes in the heap

	int null_path_length() const;
	returns the min null path length of the root node

	Type top() const;
	returns the top element of the heap
	Because this is a leftist min heap the smallest element is the root, which is the top element

	int count( Type const & ) const;
	Counts the number of object with the given paramter are currntley stored in the heap

Constructor
	Leftist_heap<Type>::Leftist_heap()
	initializes the root node to nulptr and sets the heap size to 0

Destructor
	Leftist_heap<Type>::~Leftist_heap()
	Calls the clear function on the root node which recursuvley removes elements from the heap until empty
*/
template <typename Type>
class Leftist_heap {
	private:
		Leftist_node<Type> *root_node;
		int heap_size;

	public:
		Leftist_heap();
		Leftist_heap( Leftist_heap const & );
		~Leftist_heap();

		void swap( Leftist_heap &heap );
		Leftist_heap &operator=( Leftist_heap );

		bool empty() const;
		int size() const;
		int null_path_length() const;
		Type top() const;
		int count( Type const & ) const;

		void push( Type const & );
		Type pop();
		void clear();

	// Friends

	template <typename T>
	friend std::ostream &operator<<( std::ostream &, Leftist_heap<T> const & );
};

// some sample functions are given

/*
See class definition for details
*/
template <typename Type>
Leftist_heap<Type>::Leftist_heap():
root_node( nullptr ),
heap_size( 0 ) {
	// does nothing
}

/*
Leftist_heap<Type>::Leftist_heap( Leftist_heap const &heap )
Copy constuctor for class Leftist_heap

Param 1:
	An object of type Leftist_heap which is the heap we wish to make a copy of

A pre-condition is that the functions for class Dynamic_queue are functioning correctly. 
This constuctor has no exception handling, post conditions

The constuctor copies elements of the passed in heap by using a queue to preform a bredth-first traversal of the heap
We begin by enqueing the root and then while the queue is not empty:
	dequeue the the front of the queue
	push the element which was dequeued
	enqueue the children of the node which was dequeued (we check to make sure the trees are not null)

We then deallocate memory which had been allocated by the queue
*/
template <typename Type>
Leftist_heap<Type>::Leftist_heap( Leftist_heap const &heap ):
root_node( nullptr ),
heap_size( 0 ) {
	Leftist_node<Type> *root = heap.root_node;
	Dynamic_queue< Leftist_node<Type> *> *queue = new Dynamic_queue < Leftist_node<Type> *>();
	queue->enqueue(root);
	while (!queue->empty())
	{
		Leftist_node<Type> *node = queue->dequeue();
		push(node->retrieve());
		if (node->left() != nullptr) {
			queue->enqueue(node->left()); 
		}
		if (node->right() != nullptr) {
			queue->enqueue(node->right());
		}
	}
	delete queue;
}

/*
See class definition for details
*/
template <typename Type>
Leftist_heap<Type>::~Leftist_heap() {
	clear();  // might as well use it...
}

template <typename Type>
void Leftist_heap<Type>::swap( Leftist_heap<Type> &heap ) {
	std::swap( root_node, heap.root_node );
	std::swap( heap_size, heap.heap_size );
}

template <typename Type>
Leftist_heap<Type> &Leftist_heap<Type>::operator=( Leftist_heap<Type> rhs ) {
	swap( rhs );

	return *this;
}

/*
See class definition for details
*/
template <typename Type>
bool Leftist_heap<Type>::empty() const {
	return (heap_size == 0);
}

/*
See class definition for details
*/
template <typename Type>
int Leftist_heap<Type>::size() const{
	return heap_size;
}

/*
See class definition for details
*/
template <typename Type>
int Leftist_heap<Type>::null_path_length() const {
	return root_node->null_path_length();
}

/*
See class definition for details
*/
template <typename Type>
int Leftist_heap<Type>::count(const Type &obj) const {
	return root_node->count(obj);
}

/*
See class definition for details
*/
template <typename Type>
Type Leftist_heap<Type>::top() const {
	if (empty()) {
		throw underflow();
	}
	return root_node->retrieve();
}
/*
void Leftist_heap<Type>::push (const Type &obj )
pushes a new node onto the heap

Param 1:
	An object of type Type to be pushed onto the heap

This function has no exception handling, pre or post conditions
It simply calls the push function on the root node (which is the recursive push function defined for class Leftist_node)
*/
template <typename Type>
void Leftist_heap<Type>::push (const Type &obj ) {
	Leftist_node<Type> *newNode = new Leftist_node < Type >(obj) ;
	root_node->push(newNode, root_node);
	heap_size++;
}

/*
Type Leftist_heap<Type>::pop()
pops the top elemet of the heap; because it is a leftist heap the top element is the root which is also the min element in the heap

The function has no exception handling or post condition
A pre condition for the function to work correctly is the heap must not be empty otherwise the function throws an underflow exception

We save the root node's element and then same a temp pointer to the root node
We promote the root of the left subtree to be the root of the new heap and then push the right subtree onto this new heap
finally delete memory allocated for the old root and update fields
*/
template <typename Type>
Type Leftist_heap<Type>::pop() {
	if (this == nullptr) {
		throw underflow();
	}

	Type top_node = top();
	Leftist_node<Type> *tmp = root_node;
	root_node = tmp->left();
	root_node->push(tmp->right(), root_node);
	delete tmp;
	heap_size--;

	return top_node;

}

/*
void Leftist_heap<Type>::clear()
This function clears the heap of all nodes

The is no exception handling, pre or post conditions
The function simply calls the clear function on the root node (which is the recursive clear function for class leftist_node) and updates member fields
*/
template <typename Type>
void Leftist_heap<Type>::clear() {
	root_node->clear();
	root_node = nullptr;
	heap_size = 0;
}

/*
Used to print out the heap for testing purposes
*/
template <typename Type>
std::ostream &operator<<(std::ostream &out, const Leftist_heap<Type> &heap) {
	out << "Size:     " << heap.heap_size << std::endl;

	if (!heap.empty()) {

		std::stack<Leftist_node<Type> *> traversal;

		std::stack<int> indentation;

		traversal.push(heap.root_node);
		indentation.push(0);

		while (!traversal.empty()) {
			Leftist_node<Type> *node = traversal.top();
			int indent = indentation.top();

			traversal.pop();
			indentation.pop();

			out.width(indent);
			out.fill(' ');
			out << "";

			if (node == 0) {
				out << "-" << std::endl;
			}
			else {
				out << node->retrieve() << std::endl;

				if (node->left() != 0 || node->right() != 0) {
					traversal.push(node->right());
					indentation.push(indent + 1);

					traversal.push(node->left());
					indentation.push(indent + 1);
				}
			}
		}
	}



	return out;
}

// Is an error showing up in ece250.h or elsewhere?
// Did you forget a closing '}' ?

#endif
