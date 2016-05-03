/*****************************************
 * Author: Davan Basran
 * 
 *****************************************/

#ifndef DYNAMIC_QUEUE_H
#define DYNAMIC_QUEUE_H

#ifndef nullptr
#define nullptr 0
#endif

#include <algorithm>
#include "Exception.h"
/*
Class Dynamic_queue
Memeber fields:
	int initial_capacity;
	int array_capacity;
	Type *array;
	int ihead;
	int itail;
	int entry_count;

Acessors:
	Type head() const;
	retruns the first element in the queue
	Throws underflow exception if the queue is empty

	int size() const;
	returns the number of element in the queue

	bool empty() const;
	returns true if the queue is currently empty (0 elements)

	int capacity() const;
	returns the size of the dynamically resizied array that the elements in the queue are storred in

Dynamic_queue<Type>::~Dynamic_queue()
Desructor for class Dyncamic queue
Calls delete on the array member which was dynamically allocated
*/
template <typename Type>
class Dynamic_queue {
	private:
		int initial_capacity;
		int array_capacity;
		Type *array;
		int ihead;
		int itail;
		int entry_count;

	public:
		Dynamic_queue( int = 10 );
		Dynamic_queue( Dynamic_queue const & );
		~Dynamic_queue();

		Type head() const;
		int size() const;
		bool empty() const;
		int capacity() const;

		void swap( Dynamic_queue & );
		Dynamic_queue &operator=( Dynamic_queue );
		void enqueue( Type const & );
		Type dequeue();
		void clear();

	// Friends

	template <typename T>
	friend std::ostream &operator<<( std::ostream &, Dynamic_queue<T> const & );
};

/*
Dynamic_queue<Type>::Dynamic_queue( int n)
Constuctor for class Dynamic queue

Param 1:
	integer n which is the initial capacity of the dynamic queue

There is no exception handling, pre or post conditions

Initiates:
	Allocates memory for an array of size initial_capacity
	Sets ihead, itail, entry_count to 0
*/
template <typename Type>
Dynamic_queue<Type>::Dynamic_queue( int n):
initial_capacity( std::max( n, 1 ) ),
array_capacity( initial_capacity ),
array( new Type[initial_capacity] ),
ihead( 0 ),
itail( initial_capacity - 1 ),
entry_count( 0 ) {
	itail = ihead;
}

/*
Dynamic_queue<Type>::Dynamic_queue( Dynamic_queue const &queue )
Copy Condtructor for Class Dynamic queue

Param 1:
	Object of type Dynamic_queue to be copied

There is no exception handeling, pre or post conditions

Initiates:
	ihead, itail and entry_count to that of the passed queue.
	Allocates memory for a new array to store the coppied elements

Loops through the qeueue and stores the elements in the same order and at the same indexes of the oringinal queue in
	the new copy queue
*/
template <typename Type>
Dynamic_queue<Type>::Dynamic_queue( Dynamic_queue const &queue ):
initial_capacity( queue.initial_capacity ),
array_capacity( queue.array_capacity ),
array( new Type[array_capacity] ),
ihead( queue.ihead ),
itail( queue.itail ),
entry_count( queue.entry_count ) {

	for (int i = ihead; i <= itail; i = (i +1) % capacity(), i++)
	{
		array[i] = queue.array[i];
	}
}

/*
See class definition
*/
template <typename Type>
Dynamic_queue<Type>::~Dynamic_queue() {
	delete[] array;
}

/*
See class defintion
*/
template <typename Type>
int Dynamic_queue<Type>::size() const {
	return entry_count;
}

/*
See class defintion
*/
template <typename Type>
int Dynamic_queue<Type>::capacity() const {
	return array_capacity;
}

/*
See class defintion
*/
template <typename Type>
bool Dynamic_queue<Type>::empty() const {
	return (size() == 0);
}

/*
See class defintion
*/
template <typename Type>
Type Dynamic_queue<Type>::head() const {
	if (empty()) {
		throw underflow();
	}
	return array[ihead];
}

template <typename Type>
void Dynamic_queue<Type>::swap( Dynamic_queue<Type> &queue ) {
	std::swap( initial_capacity, queue.initial_capacity );
	std::swap( array_capacity, queue.array_capacity );
	std::swap( array, queue.array );
	std::swap( ihead, queue.ihead );
	std::swap( itail, queue.itail );
	std::swap( entry_count, queue.entry_count );
}

template <typename Type>
Dynamic_queue<Type> &Dynamic_queue<Type>::operator=( Dynamic_queue<Type> rhs ) {
	swap( rhs );
	
	return *this;
}

/*
void Dynamic_queue<Type>::enqueue( Type const &obj )
Pushes an object to the back of the queue

Param 1:
	an Object of type Type to be enqueued

There is no exception handling, pre or post conditions

2 cases:
	Case 1:
	If we have filled the array, then double the size of the array
	proceed to copy entries over in a efficient way (index 0 to n) and deallocate memory used by old array and push element

	Case 2:
	update the index of itail by using a mod counter then insert the object in the new position
*/
template <typename Type>
void Dynamic_queue<Type>::enqueue( Type const &obj ) {
	if ((itail == capacity() - 1) && (size()!= 0)) {
		// double the capacity
		Type *newArray = new Type[capacity() * 2];
		int j = ihead;
		for (int i = 0; i < size(); j = (j + 1) % capacity(), i++)
		{
			newArray[j] = array[i];
		}
		delete[] array;
		array = newArray;
		
		array_capacity = capacity() * 2;

		// update tail and head for new array size
		ihead = 0;
		itail = (ihead + size() ) % capacity();
		array[itail] = obj;
	}
	else {
		itail = (ihead + size()) % capacity();
		array[itail] = obj;
	}
	entry_count++;
}

/*
Type Dynamic_queue<Type>::dequeue()
Pops an element from the font of the queue

There is no exception handling, pre or post conditions

2 cases:
	Case 1:
	After popping the element, the number of elements in the array is not 1/4 of the max capacity and the capacity is not under 20

	Case 2:
	After popping the element the number of elements of 1/4 of the capacity and the capacity is over 19
	Allocate memory of an array that is half the size of original array and proceed to copy elements over efficently (index 0 to n)
	deallocate memory of old array
	update relevent members
*/
template <typename Type>
Type Dynamic_queue<Type>::dequeue() {
	if ( empty() ) {
		throw underflow();
	}

	Type obj = head();
	entry_count--;
	ihead = (ihead +1) % capacity() ;

	if ( (size() == (capacity() / 4) ) && 
		( (capacity() / 2) >= initial_capacity) ) {
		Type *newArray = new Type[capacity() / 2];
		//copy entires from old array over new entires will begin at index 0;
		int j = 0;
		for (int i = ihead; i <= itail; i = (i + 1) % capacity(), j++)
		{
			newArray[j] = array[i];
		}
		delete[] array;
		array = newArray;
		//update size for new array
		array_capacity = capacity() / 2;

		// update tail and head for new array size
		ihead = 0;
		itail = (ihead + size() - 1) % capacity();
	}
	return obj;
}

/*
void Dynamic_queue<Type>::clear()
Deletes old array and replaces it by a new empty array

Assumes memory for array has been allocated before calling function
No exception handling, or post conditions

Allocate memory of new empty array and deallocate memory used by old array
*/
template <typename Type>
void Dynamic_queue<Type>::clear() {
	Type *newArray = new Type[initial_capacity];
	ihead = 0;
	itail = initial_capacity - 1;
	entry_count = 0;
	delete[] array;
	array = newArray;
}

// You can modify this function however you want:  it will not be tested

template <typename Type>
std::ostream &operator<<( std::ostream &out, Dynamic_queue<Type> const &queue ) {
	// I don't know how you are implementing your queue so I cannot print it.
	// If you want, you can print whatever you want here and then call cout
	// in the driver.

	// Remember to redirect to out, e.g.,
	//      out << "Hello!";

	return out;
}

// Is an error showing up in ece250.h or elsewhere?
// Did you forget a closing '}' ?

#endif
