
////////////////////////
// print_collection.h //
////////////////////////

#ifndef print_collection_h
#define print_collection_h

#include <ostream>

template <typename TCollection>
std::ostream & print_collection(std::ostream & out, const TCollection & collection)
{
    bool first = true;

    out << '{';

    for (const typename TCollection::value_type & element : collection)
    {
        if (first)
        {
            first = false;
        }
        else
        {
            out << ", ";
        }

        out << element;
    }

    out << '}';

    return out;
}

#endif // print_collection_h
