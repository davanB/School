/*****************************************
 *  Author: Davan Basran
 *****************************************/

#ifndef DOUBLE_HASH_TABLE_H
#define DOUBLE_HASH_TABLE_H

#ifndef nullptr
#define nullptr 0
#endif

#include "Exception.h"

/*
enum bin_state_t
used to keep track of the current state of a bin in the hash table
3 states:
	UNOCCUPIED //there is no element in the bin and uppon accessed will stop bin searching
	OCCUPIED // there is currently an element saved in the bin, when accessed will continue bin searching
	ERASED // an occuiped element has been lazy erased and can no longer be accessed, will continue when accesed on bin search
*/
enum bin_state_t { UNOCCUPIED, OCCUPIED, ERASED };

/*
class Quadratic_hash_table
Stores elements in a hashtable to be accessed in ideally order 1 time

member fields
	int count; // number of elements in hash table (occupied bins)
	int power; // represent 2^power which controls the inital size of the hash table (a power of 2)
	int array_size; // is he size of the hash table
	int mask; // The number of bits we want to keep when using logical AND for modulo operations (bitwise)
	double erased_bins; // added count to keep track of erased bins
	Type *array; // an array pointer which is the hashtable in memory
	bin_state_t *occupied; //an array which keeps track of the state of the bins

Accessors and Functions:
	int hash( Type const & ) const;
	Functions takes in a single parameter which is the object we wish to hash and map onto a bin in the hash table

	void print() const;
	itterates through the hash table and prints the objects stored in it

	int size() const;
	returns the number of elements currently stored in the hash table (are occupied bins)

	int capacity() const;
	returns the size of the array member field which represents the size of the hash table

	double load_factor() const;
	calculates the load factor of the hash table (object occupied + erased)/total bins

	bool empty() const;
	returns true if the hash table is empty

	bool member( Type const & ) const;
	returns true of the passed object is in the hash table (in an occupied bin)

	Type bin( int ) const;
	returns the object stored at the passed in bin number

Constructor:
	Quadratic_hash_table<Type>::Quadratic_hash_table( int m )
	Takes as a paramter the power member field which dictates the intial size of the array (2^power)
	initializes the member fields to the defaults and allocares memory for the 2 arrays
	Initializes all bins to the UNOCCUPIED states

Destructor:
	~Quadratic_hash_table();
	deallocates memory for the hash table stored as an array and the array which keeps track of the bin states
*/
template <typename Type>
class Quadratic_hash_table {
	private:
		int count;
		int power;
		int array_size;
		int mask;
		double erased_bins; 
		Type *array;
		bin_state_t *occupied;

		int hash( Type const & ) const;

	public:
		Quadratic_hash_table( int = 5 );
		~Quadratic_hash_table();
		int size() const;
		int capacity() const;
		double load_factor() const;
		bool empty() const;
		bool member( Type const & ) const;
		Type bin( int ) const;

		void print() const;

		void insert( Type const & );
		bool erase( Type const & );
		void clear();

	// Friends

	template <typename T>
	friend std::ostream &operator<<( std::ostream &, Quadratic_hash_table<T> const & );
};

/*
See class definition for details
*/
template <typename Type>
Quadratic_hash_table<Type>::Quadratic_hash_table( int m ):
count( 0 ), power( m ),
array_size( 1 << power ),
mask( array_size - 1 ),
array( new Type[array_size] ),
occupied( new bin_state_t[array_size] ) {
	erased_bins = 0;
	for ( int i = 0; i < array_size; ++i ) {
		occupied[i] = UNOCCUPIED;
	}
}

/*
See class definition for details
*/
template <typename Type>
Quadratic_hash_table<Type>::~Quadratic_hash_table() {
	delete[] array;
	delete[] occupied;
}

/*
See class definition for details
*/
template <typename Type>
int Quadratic_hash_table<Type>::size() const {
	return count;
}

/*
See class definition for details
*/
template <typename Type> 
int Quadratic_hash_table<Type>::capacity() const {
	return array_size;
}

/*
See class definition for details
*/
template <typename Type>
double Quadratic_hash_table<Type>::load_factor() const {
	return (static_cast<double>(size()) + erased_bins) /
		static_cast<double>(capacity());
}

/*
See class definition for details
*/
template <typename Type>
bool Quadratic_hash_table<Type>::empty() const {
	return (size() == 0);
}

/*
See class definition for details
*/
template <typename Type>
bool Quadratic_hash_table<Type>::member(Type const &obj) const {
	int bin_number = hash(obj);
	if (occupied[bin_number] == UNOCCUPIED) return false;
	for (int k = 0; k < capacity(); k++) 
	{
		if (array[bin_number] == obj) {
			if (occupied[bin_number] == OCCUPIED)
				return true;
			else
				return false;
		}
		bin_number = (bin_number + k) & mask;
	}
	//iterated through whole array
	return false;
}

/*
See class definition for details
*/
template <typename Type>
Type Quadratic_hash_table<Type>::bin(int n) const {
	return array[n];
}

/*
See class definition for details
*/
template <typename Type>
void Quadratic_hash_table<Type>::print() const {
	for (int i = 0; i < capacity(); i++)
	{
		if (occupied[i] == OCCUPIED) {
			std::cout << "[" <<  array[i] << "]";
		}
		else if (occupied[i] == ERASED) {
			std::cout << "[X]";
		}
		std::cout << '[ ]';
	}
}

/*
See class definition for details
*/
template <typename Type>
int Quadratic_hash_table<Type>::hash(Type const &obj) const {
	int bin_number = static_cast<int>(obj) & mask;
	return (bin_number < 0 ? bin_number + capacity() : bin_number);
}

/*
void Quadratic_hash_table<Type>::insert(Type const &obj)
Inserts a new object into the hash table

Param1:
	An object of type Type to be stored in the hash table

Precondition: the hash function is working correctly and maps the obj to the correct bin
this function does not have any post conditions and no exception handling

We first map the object to a bin
We then inset the obj into the bin
	If the bin is occupied we use Quadratic probing to find a new bin to store the obj in
*/
template <typename Type>
void Quadratic_hash_table<Type>::insert(Type const &obj) {
	if (size() == capacity()) throw overflow();
	if (!member(obj)) { //check if obj is in hash if not add it in
		int bin_number = hash(obj);
		for (int i = 0; i < capacity(); i++)
		{
			if (occupied[bin_number] != OCCUPIED) {
				if (occupied[bin_number] == ERASED) {
					erased_bins--;
				}
				array[bin_number] = obj;
				occupied[bin_number] = OCCUPIED;
				count++;
				return;
			}
			else {
				bin_number = (bin_number + i) & mask;
			}
		}
	}
}

/*
bool Quadratic_hash_table<Type>::erase(Type const &obj)
uses Lazy deletion to erase an object from the hash table, returns true if th obj has been erased

Param1:
	An object of type Type to be erased from the hash table

Preconditions:
	The member funciton works correctly so that it is not assumed passed obj is in the hash table or not
	The hash function works correctly
There are no post conditions and does not have exception handling

Searching the bin which would contain the obj and sets the bin to the ERASED state and updates appropriate variables
If the obj is not in the bin we use Quadratic probing to find the obj 
*/
template <typename Type>
bool Quadratic_hash_table<Type>::erase(Type const &obj) {
	if (!member(obj)) {
		return false;
	}
	else {
		int bin_number = hash(obj);
		for (int i = 0; i < capacity(); i++)
		{
			if (array[bin_number] == obj) {
				occupied[bin_number] = ERASED;
				count--;
				erased_bins++;
				return true;
			}
			else {
				bin_number = (bin_number + i) & mask;
			}
		}
		return false;
	}
}

/*
void Quadratic_hash_table<Type>::clear()
resets the hash table

There are no pre, post conditions and no exception handling
Sets all bins to UNOCCUPIED and resets variables to the defaults
*/
template <typename Type>
void Quadratic_hash_table<Type>::clear() {
	for (int i = 0; i < array_size; ++i) {
		occupied[i] = UNOCCUPIED;
	}
	erased_bins = 0;
	count = 0;
}

template <typename T>
std::ostream &operator<<( std::ostream &out, Quadratic_hash_table<T> const &hash ) {
	for ( int i = 0; i < hash.capacity(); ++i ) {
		if ( hash.occupied[i] == UNOCCUPIED ) {
			out << "- ";
		} else if ( hash.occupied[i] == ERASED ) {
			out << "x ";
		} else {
			out << hash.array[i] << ' ';
		}
	}

	return out;
}

#endif
