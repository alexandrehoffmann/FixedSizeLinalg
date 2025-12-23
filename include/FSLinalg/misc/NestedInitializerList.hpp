#ifndef FSLINALG_MISC_NESTED_INITIALIZER_LIST_HPP
#define FSLINALG_MISC_NESTED_INITIALIZER_LIST_HPP

namespace FSLinalg
{
namespace misc
{

template<typename T, size_t rank>
struct NestedInitializerListTraits
{
    static_assert(rank > 0);

    using Type = std::initializer_list<typename NestedInitializerListTraits<T, rank-1>::Type>;
};

template<typename T>
struct NestedInitializerListTraits<T,1>
{
    using Type = std::initializer_list<T>;
};

template<typename T, size_t rank>
using NestedInitializerList = typename NestedInitializerListTraits<T,rank>::Type;

} //namespace misc
} // namespace FSLinalg

#endif // FSLINALG_MISC_NESTED_INITIALIZER_LIST_HPP
