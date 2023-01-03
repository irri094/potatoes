#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
template<typename T>
using ordered_set = tree<T,null_type,less<T>,
rb_tree_tag, tree_order_statistics_node_update>;
//find_by_order(k):iterator to kth smallest(0 ind)
//order_of_key(x) :  no of items  < x
