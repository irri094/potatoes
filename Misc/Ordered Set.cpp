#include <ext/pb_ds/assoc_container.hpp>
#include <ext/pb_ds/tree_policy.hpp>
using namespace __gnu_pbds;
typedef tree<int,null_type,less<int>,rb_tree_tag,
tree_order_statistics_node_update > ordered_set;
//find_by_order(k):iterator to kth smallest(0 ind)
//order_of_key(x) :  no of items  <= x
